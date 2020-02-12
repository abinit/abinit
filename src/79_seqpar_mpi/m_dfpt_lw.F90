!!****m* ABINIT/m_dfpt_lw
!! NAME
!!  m_dfpt_lw
!!
!! FUNCTION
!!  Calculation of spatial dispertion magnitudes (quadrupole and
!!  flexoelectric tensor) from the DFPT long-wave approach.
!!
!! COPYRIGHT
!!  Copyright (C) 2018 ABINIT group (MR,MS)
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

module m_dfpt_lw
    
 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_dtset
 use m_dtfil
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_wfk
 use m_nctk
 use m_hamiltonian
 use m_ebands
 use m_pawrhoij
 use m_pawtab
 use m_wffile

 use m_hdr
 use m_io_tools,   only : file_exists
 use m_ioarr,      only : read_rhor, fftdatar_write_from_hdr
 use m_wfk
 use m_rf2,        only : rf2_getidirs
 use m_kg,         only : getcut, getph, getmpw, kpgio
 use m_abicore,    only : appdig
 use m_fft,        only : fourdp
 use m_dfpt_rhotov, only : dfpt_rhotov
 use m_spacepar,   only : hartredq, setsym
 use m_cgtools,    only : dotprod_vn
 use m_symkpt,     only : symkpt
 use m_mpinfo,     only : distrb2, initmpi_band, proc_distrb_cycle
 use m_initylmg,   only : initylmg
 use m_inwffil,    only : inwffil
 use m_dfpt_lwwf
 use m_dynmat,     only : cart39

 implicit none

 private
!***

 public :: dfpt_qdrpole
 public :: dfpt_flexo

! *************************************************************************

contains 

!!****f* ABINIT/dfpt_qdrpole
!! NAME
!!  dfpt_qdrpole
!!
!! FUNCTION
!!  This routine computes the elements of the quadrupole and the P^(1) tensor as the 
!!  second order q-gradient (d2/dq1dq2) of the charge response to an atomic 
!!  displacement (see, e.g., M.Royo and M. Stengel to be published).
!!  The atoms and atomic displacements directions that are evaluated 
!!  are fixed by the rfatpol and rfdir variables in the input, rfdir also
!!  fixes the directions for the dq2 derivative, whereas dq1 is evaluated
!!  in all three directions.
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2017 ABINIT group (MR,MS)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  codvsn=code version
!!  doccde(mband*nkpt*nsppol)=derivative of occupancies wrt the energy
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  kxc(nfft,nkxc)=exchange and correlation kernel
!!  mpert=maximum number of perturbations for output processing
!!  mpi_enreg=information about MPI parallelization
!!  nattyp(ntypat)= # atoms of each type.
!!  nfft=(effective) number of FFT grid points (for this proc)
!!  ngfft(1:18)=integer array with FFT box dimensions and other
!!  nkpt=number of k points in the full BZ
!!  nkxc=second dimension of the kxc array. If /=0, the XC kernel must be computed.
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(mband*nkpt*nsppol)=occup number for each band (often 2) at each k point
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data for the GS (DUMMY)
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data (DUMMY)
!!  pertsy(3,natom+6)=set of perturbations that form a basis for all other perturbations
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  rhog(2,nfftf)=array for Fourier transform of GS electron density
!!  rhor(nfftf,nspden)=array for GS electron density in electrons/bohr**3.
!!  timrev=1 if time-reversal preserves the q wavevector; 0 otherwise.
!!  ucvol=unit cell volume in bohr**3.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  blkflg(3,mpert,3,mpert,3,mpert)= ( 1 if the element of the 3dte
!!   has been calculated ; 0 otherwise )
!!  d3etot(2,3,mpert,3,mpert,3,mpert)= matrix of the 3DTE
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  For historical reasons the Quadrupole and P^(1) tensor has been obtained from
!!  the q-gradient of E^{\tau_{\kappa\beta}* \epsilon_{alpha}}. This is 
!!  different from the equations appearing in PRX 9, 021050 (2019) which
!!  stand for the q-gradient of E^{\epsilon_{alpha}* \tau_{\kappa\beta}}.
!!  By this reason there is an additional -1 factor here when computing 
!!  the physical observables. 
!! 
!! PARENTS
!!   
!!  respfn
!!
!! CHILDREN
!!
!!  appdig, distrb2, dfpt_qdrpout, dfpt_qdrpwf, dfpt_rhotov, distrb2, dotprod_vn, ebands_init, 
!!  ebands_free, fftdatar_write_from_hdr, 
!!  fourdp, getcut, getmpw, getph, hdr_init, hdr_update, init_hamiltonian, 
!!  initmpi_band, initylmg, inwffil, kpgio, load_spin_hamiltonian,
!!  hartredq, read_rhor, rf2_getidirs, setsym, symkpt, WffClose,  
!!  wfk_open_read, wrtout, xmpi_sum
!!
!! SOURCE

subroutine dfpt_qdrpole(atindx,blkflg,codvsn,d3etot,doccde,dtfil,dtset,&
&          gmet,gprimd,kxc,mpert,&
&          mpi_enreg,nattyp,nfft,ngfft,nkpt,nkxc,&
&          nspden,nsppol,occ,pawrhoij,pawtab,pertsy,psps,rmet,rprimd,rhog,rhor,&
&          timrev,ucvol,xred)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_qdrpole'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpert,nfft,nkpt,nkxc,nspden,nsppol,timrev
 integer,intent(inout) :: blkflg(3,mpert,3,mpert,3,mpert)
 real(dp),intent(inout) :: d3etot(2,3,mpert,3,mpert,3,mpert)
 real(dp),intent(in) :: ucvol
 character(len=8), intent(in) :: codvsn
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(hdr_type) :: hdr_den
 type(pseudopotential_type),intent(in) :: psps
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)
!arrays
 integer,intent(in) :: atindx(dtset%natom)
 integer,intent(in) :: nattyp(dtset%ntypat),ngfft(18)
 integer,intent(in) :: pertsy(3,dtset%natom+6)
 real(dp),intent(in) :: doccde(dtset%mband*nkpt*dtset%nsppol)
 real(dp),intent(in) :: gmet(3,3), gprimd(3,3)
 real(dp),intent(in) :: kxc(nfft,nkxc)
 real(dp),intent(in) :: occ(dtset%mband*nkpt*dtset%nsppol)
 real(dp),intent(in) :: rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: rhog(2,nfft),rhor(nfft,nspden)
 real(dp),intent(in) :: xred(3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer :: ask_accurate,bdtot_index,bdtot1_index,bantot_rbz
 integer :: cplex,formeig,forunit,gscase,iatpert,iatpert_cnt,iatpol
 integer :: iatdir,icg,ierr,ii,ikg,ikpt,ikpt1,ilm,iq1dir,iq1grad,iq1grad_cnt
 integer :: iq1q2grad,iq1q2grad_var,iq2dir,iq2grad,iq2grad_cnt,ireadwf0,isppol,istwf_k
 integer :: jj,master,matom,matpert,mcg,me,mgfft
 integer :: mkmem_rbz,mpw,my_nkpt_rbz
 integer :: natpert,nband_k,nfftot,nhat1grdim,nkpt_rbz,npw_k,npw1_k
 integer :: nq1grad,nq1q2grad,nq2grad,nsym1,nylmgr,n3xccc
 integer :: optorth,optene,option,optres
 integer :: pawread,pertcase,qdir,spaceworld
 integer :: usexcnhat,useylmgr
 integer,parameter :: formeig1=1
 integer,parameter :: re=1,im=2
 real(dp) :: boxcut,doti,dotr,dum_scl,ecut_eff,ecut,etotal,fermie,gsqcut,residm
 real(dp) :: vres2, wtk_k
 logical :: non_magnetic_xc
 character(len=500) :: msg                   
 character(len=fnlen) :: fi1o,fiwfatdis,fiwfefield,fiwfddk,fiwfdkdk
 type(ebands_t) :: bs_rbz
 type(hdr_type) :: hdr0
 type(wvl_data) :: wvl 
 type(wffile_type) :: wffgs,wfftgs
 type(gs_hamiltonian_type) :: gs_hamkq

!arrays
 integer,allocatable :: bz2ibz_smap(:,:)
 integer,allocatable :: indkpt1(:), indkpt1_tmp(:),indsy1(:,:,:),irrzon1(:,:,:)
 integer,allocatable :: istwfk_rbz(:),kg(:,:),kg_k(:,:)
 integer,allocatable :: nband_rbz(:),npwarr(:),npwtot(:)
 integer,allocatable :: pert_atdis(:,:), pert_atdis_tmp(:,:)
 integer,allocatable :: q1grad(:,:),q1grad_tmp(:,:),q1q2grad(:,:),q2grad(:,:),q2grad_tmp(:,:)
 integer,allocatable :: qdrflg(:,:,:,:)
 integer,allocatable :: symaf1(:),symrc1(:,:,:),symrl1(:,:,:)
 real(dp) :: kpoint(3)
 real(dp),allocatable :: cg(:,:),doccde_rbz(:)
 real(dp),allocatable :: eigen0(:),eqgradhart(:,:,:,:)
 real(dp),allocatable :: kpt_rbz(:,:)
 real(dp),allocatable :: nhat(:,:),nhat1(:,:),nhat1gr(:,:,:) 
 real(dp),allocatable :: occ_k(:),occ_rbz(:)
 real(dp),allocatable :: ph1d(:,:),phnons1(:,:,:) 
 real(dp),allocatable :: qdrpwf(:,:,:,:),qdrpwf_k(:,:,:,:)
 real(dp),allocatable :: qdrpwf_t1(:,:,:,:),qdrpwf_t1_k(:,:,:,:)
 real(dp),allocatable :: qdrpwf_t2(:,:,:,:),qdrpwf_t2_k(:,:,:,:)
 real(dp),allocatable :: qdrpwf_t3(:,:,:,:),qdrpwf_t3_k(:,:,:,:)
 real(dp),allocatable :: qdrpwf_t4(:,:,:,:),qdrpwf_t4_k(:,:,:,:)
 real(dp),allocatable :: qdrpwf_t5(:,:,:,:),qdrpwf_t5_k(:,:,:,:)
 real(dp),allocatable :: rhog1_atdis(:,:,:)
 real(dp),allocatable :: rhog1_tmp(:,:)
 real(dp),allocatable :: rhor1_efield(:,:,:),rhor1_tmp(:,:),rhor1_real(:,:)
 real(dp),allocatable :: tnons1(:,:)
 real(dp),allocatable :: vhartr1(:),vhxc1_atdis(:,:),vhxc1_efield(:,:)
 real(dp),allocatable :: vpsp1(:),vqgradhart(:),vresid1(:,:),vxc1(:,:)
 real(dp),allocatable :: dum_vxc(:,:)
 real(dp),allocatable :: wtk_folded(:), wtk_rbz(:)
 real(dp),allocatable,target :: vtrial1(:,:)
 real(dp),allocatable :: xccc3d1(:)
 real(dp),allocatable :: ylm(:,:),ylm_k(:,:),ylmgr(:,:,:),ylmgr_k(:,:,:)
 type(pawrhoij_type),allocatable :: pawrhoij_read(:)
 type(wfk_t),allocatable :: wfk_t_atdis(:),wfk_t_efield(:),wfk_t_ddk(:),wfk_t_dkdk(:)
 
! *************************************************************************

 DBG_ENTER("COLL")

!Anounce start of quadrupole tensor calculation
 write(msg, '(a,80a,a,a,a)' ) ch10,('=',ii=1,80),ch10,&
&   ' ==> Compute Quadrupole Tensor <== ',ch10
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

!Keep track of total time spent in dfpt_qdrpole
!!To be implemented if required

!Init parallelism
 spaceworld=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt
 master=0

!Get FFT grid(s) sizes (be careful !) See NOTES in the comments at the beginning of respfn.F90
 mgfft=dtset%mgfft
 ecut=dtset%ecut
 ecut_eff=ecut*(dtset%dilatmx)**2

!Compute large sphere cut-off gsqcut
 call getcut(boxcut,ecut,gmet,gsqcut,dtset%iboxcut,std_out,dtset%qptn,ngfft)

!Various initializations
 cplex=2-timrev
 matom=dtset%natom
 usexcnhat=0
 n3xccc=0
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 non_magnetic_xc=.true.
 

!Generate the 1-dimensional phases
 ABI_ALLOCATE(ph1d,(2,3*(2*dtset%mgfft+1)*dtset%natom))
 call getph(atindx,dtset%natom,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),ph1d,xred)

!################# PERTURBATIONS AND q-GRADIENTS LABELLING ############################

!Determine which atomic displacements and q-gradient directions have to be evaluated
!taking into account the perturbation symmetries
 matpert=dtset%natom*3
 ABI_ALLOCATE(pert_atdis_tmp,(3,matpert))
 pert_atdis_tmp=0
 iatpert_cnt=0
 do iatpol=1,matom 
   do iatdir=1,3
     if (pertsy(iatdir,iatpol)==1) then
       iatpert_cnt=iatpert_cnt+1
       pert_atdis_tmp(1,iatpert_cnt)=iatpol                !atom displaced
       pert_atdis_tmp(2,iatpert_cnt)=iatdir                !direction of displacement
       pert_atdis_tmp(3,iatpert_cnt)=(iatpol-1)*3+iatdir   !like pertcase in dfpt_loopert.f90
     end if
   end do
 end do
 natpert=iatpert_cnt
 ABI_ALLOCATE(pert_atdis,(3,natpert))
 do iatpert=1,natpert
   pert_atdis(:,iatpert)=pert_atdis_tmp(:,iatpert)
 end do
 ABI_DEALLOCATE(pert_atdis_tmp)
 
 !The q2grad is related with the response to the electric field
 ABI_ALLOCATE(q2grad_tmp,(3,3))
 q2grad_tmp=0
 iq2grad_cnt=0
 do iq2grad=1,3
   if (pertsy(iq2grad,dtset%natom+2)==1) then
     iq2grad_cnt=iq2grad_cnt+1
     q2grad_tmp(1,iq2grad_cnt)=dtset%natom+2               !electric field pert
     q2grad_tmp(2,iq2grad_cnt)=iq2grad                     !electric field direction
     q2grad_tmp(3,iq2grad_cnt)=matpert+3+iq2grad           !like pertcase in dfpt_loopert.f90
   end if
 end do
 nq2grad=iq2grad_cnt
 ABI_ALLOCATE(q2grad,(3,nq2grad))
 do iq2grad=1,nq2grad
   q2grad(:,iq2grad)=q2grad_tmp(:,iq2grad)
 end do
 ABI_DEALLOCATE(q2grad_tmp)

 !The q1grad is related with the response to the ddk 
 ABI_ALLOCATE(q1grad_tmp,(3,3))
 q1grad_tmp=0
 iq1grad_cnt=0
 do iq1grad=1,3
   if (pertsy(iq1grad,dtset%natom+1)==1) then
     iq1grad_cnt=iq1grad_cnt+1
     q1grad_tmp(1,iq1grad_cnt)=dtset%natom+1                         !ddk perturbation
     q1grad_tmp(2,iq1grad_cnt)=iq1grad                               !ddk or ddq direction
     q1grad_tmp(3,iq1grad_cnt)=matpert+iq1grad                       !like pertcase in dfpt_loopert.f90
   end if
 end do
 nq1grad=iq1grad_cnt
 ABI_ALLOCATE(q1grad,(3,nq1grad))
 do iq1grad=1,nq1grad
   q1grad(:,iq1grad)=q1grad_tmp(:,iq1grad)
 end do
 ABI_DEALLOCATE(q1grad_tmp)

 !For the evaluation of the 2nd order q-gradient, the 9 directios are activated because 
 !currently the program calculates by defect all the components of the d2_dkdk perturbation.
 !TODO: This will have to be modified in the future when ABINIT enables to calculate specific
 !components of the d2_dkdk
 nq1q2grad=9
 ABI_ALLOCATE(q1q2grad,(4,nq1q2grad))
 do iq1q2grad=1,nq1q2grad
   call rf2_getidirs(iq1q2grad,iq1dir,iq2dir)
   if (iq1dir==iq2dir) then
     q1q2grad(1,iq1q2grad)=dtset%natom+10                    !d2_dkdk perturbation diagonal elements
   else
     q1q2grad(1,iq1q2grad)=dtset%natom+11                    !d2_dkdk perturbation off-diagonal elements
   end if
   q1q2grad(2,iq1q2grad)=iq1dir                              !dq1 direction 
   q1q2grad(3,iq1q2grad)=iq2dir                              !dq2 direction 
   iq1q2grad_var=iq1q2grad
   if (iq1q2grad>6) iq1q2grad_var=iq1q2grad-3                !Lower=Upper diagonal triangle matrix
   q1q2grad(4,iq1q2grad)=iq1q2grad_var+(dtset%natom+6)*3     !like pertcase in dfpt_loopert.f90
 end do

!################# ELECTROSTATIC CONTRIBUTIONS  #######################################

!This is necessary to deactivate paw options in the dfpt_rhotov routine
 ABI_DATATYPE_ALLOCATE(pawrhoij_read,(0))
 pawread=0
 nhat1grdim=0
 ABI_ALLOCATE(nhat1gr,(0,0,0))
 ABI_ALLOCATE(nhat,(nfft,nspden))
 nhat=zero
 ABI_ALLOCATE(nhat1,(cplex*nfft,nspden))
 nhat1=zero

!Read the electric field density response from a disk file(rhor1_efield), calculates the FFT 
!(rhog1_tmp) and the first order Hartree and xc potentials(vhxc1_efield). 
!TODO: In the call to read_rhor there is a security option that compares with the header
!hdr. Not activated at this moment.
 ABI_ALLOCATE(rhog1_tmp,(2,nfft))
 ABI_ALLOCATE(rhor1_efield,(nq2grad,cplex*nfft,nspden))
 ABI_ALLOCATE(rhor1_tmp,(cplex*nfft,nspden))
 ABI_ALLOCATE(rhor1_real,(1*nfft,nspden))
 ABI_ALLOCATE(vhartr1,(cplex*nfft))
 ABI_ALLOCATE(vhxc1_efield,(nq2grad,cplex*nfft))
 ABI_ALLOCATE(vpsp1,(cplex*nfft))
 ABI_ALLOCATE(vtrial1,(cplex*nfft,nspden))
 ABI_ALLOCATE(vresid1,(cplex*nfft,nspden))
 ABI_ALLOCATE(vxc1,(cplex*nfft,nspden))
 ABI_ALLOCATE(dum_vxc,(nfft,nspden))
 ABI_ALLOCATE(xccc3d1,(cplex*n3xccc))
 vpsp1=zero; vtrial1=zero; dum_vxc=zero
 optene=0; optres=1 
 do iq2grad=1,nq2grad
   pertcase=q2grad(3,iq2grad)
   
   !Reads a real first order density
   call appdig(pertcase,dtfil%fildens1in,fi1o)
   call read_rhor(fi1o, 1, nspden, nfft, ngfft, pawread, mpi_enreg, rhor1_real, &
    & hdr_den, pawrhoij_read, spaceworld)

   !Perform FFT rhor1 to rhog1
   call fourdp(cplex,rhog1_tmp,rhor1_real,-1,mpi_enreg,nfft,1,ngfft,0)

   !Accumulate density in meaningful complex arrays
   if (timrev==0) then
     do ii=1,nfft
       jj=ii*2
       rhor1_tmp(jj-1,:)=rhor1_real(ii,:)
     end do
   else if (timrev==1) then
     rhor1_tmp(:,:)=rhor1_real(:,:)
   end if
   rhor1_efield(iq2grad,:,:)=rhor1_tmp(:,:)

   !Calculate first order Hartree and xc potentials
   call dfpt_rhotov(cplex,dum_scl,dum_scl,dum_scl,dum_scl,dum_scl, &
    & gsqcut,q2grad(2,iq2grad),dtset%natom+2,&
    & dtset%ixc,kxc,mpi_enreg,dtset%natom,nfft,ngfft,nhat,nhat1,nhat1gr,nhat1grdim,nkxc,&
    & nspden,n3xccc,non_magnetic_xc,optene,optres,dtset%qptn,rhog,rhog1_tmp,rhor,rhor1_tmp,&
    & rprimd,ucvol,psps%usepaw,usexcnhat,vhartr1,vpsp1,vresid1,vres2,vtrial1,dum_vxc,vxc1,&
    & xccc3d1,dtset%ixcrot)


   !Accumulate the potential in meaningful arrays
   vhxc1_efield(iq2grad,:)=vtrial1(:,nspden)

 end do

!Read the atomic displacement density response from a disk file, calculate the FFT 
!(rhog1_atdis) and the first order Hartree and xc potentials(vhxc1_atdis). 
 ABI_ALLOCATE(rhog1_atdis,(natpert,2,nfft))
 ABI_ALLOCATE(vhxc1_atdis,(natpert,cplex*nfft))
 vtrial1=zero
 do iatpert= 1, natpert
   iatpol=pert_atdis(1,iatpert)
   iatdir=pert_atdis(2,iatpert)
   pertcase=pert_atdis(3,iatpert)

   !Reads a real first order density
   call appdig(pertcase,dtfil%fildens1in,fi1o)
   call read_rhor(fi1o, 1, nspden, nfft, ngfft, pawread, mpi_enreg, rhor1_real, &
    & hdr_den, pawrhoij_read, spaceworld)

   !Perform FFT rhor1 to rhog1
   call fourdp(cplex,rhog1_tmp,rhor1_real,-1,mpi_enreg,nfft,1,ngfft,0)

   !Accumulate density in meaningful complex arrays
   if (timrev==0) then
     do ii=1,nfft
       jj=ii*2
       rhor1_tmp(jj-1,:)=rhor1_real(ii,:)
     end do
   else if (timrev==1) then
     rhor1_tmp(:,:)=rhor1_real(:,:)
   end if
   rhog1_atdis(iatpert,:,:)=rhog1_tmp(:,:)

   !Calculate first order Hartree and xc potentials
   call dfpt_rhotov(cplex,dum_scl,dum_scl,dum_scl,dum_scl,dum_scl, &
    & gsqcut,iatdir,iatpol,&
    & dtset%ixc,kxc,mpi_enreg,dtset%natom,nfft,ngfft,nhat,nhat1,nhat1gr,nhat1grdim,nkxc,&
    & nspden,n3xccc,non_magnetic_xc,optene,optres,dtset%qptn,rhog,rhog1_tmp,rhor,rhor1_tmp,&
    & rprimd,ucvol,psps%usepaw,usexcnhat,vhartr1,vpsp1,vresid1,vres2,vtrial1,dum_vxc,vxc1,&
    & xccc3d1,dtset%ixcrot)

   !Accumulate the potential in meaningful arrays
   vhxc1_atdis(iatpert,:)=vtrial1(:,nspden)

 end do

 !These arrays will not be used anymore (for the moment)
 ABI_DEALLOCATE(rhor1_real)
 ABI_DEALLOCATE(rhor1_tmp)
 ABI_DEALLOCATE(vhartr1)
 ABI_DEALLOCATE(vpsp1)
 ABI_DEALLOCATE(vtrial1)
 ABI_DEALLOCATE(vresid1)
 ABI_DEALLOCATE(vxc1)
 ABI_DEALLOCATE(dum_vxc)
 ABI_DEALLOCATE(xccc3d1)

 ABI_DATATYPE_DEALLOCATE(pawrhoij_read)
 ABI_DEALLOCATE(nhat1gr)
 ABI_DEALLOCATE(nhat)
 ABI_DEALLOCATE(nhat1)

!Calculate the electrostatic contribution from the q-gradient of the Hartree potential
 ABI_ALLOCATE(vqgradhart,(2*nfft))
 ABI_ALLOCATE(rhor1_tmp,(2*nfft,nspden))
 ABI_ALLOCATE(eqgradhart,(2,natpert,nq2grad,nq1grad))
 ABI_ALLOCATE(qdrflg,(matom,3,3,3))
 qdrflg=0
 rhor1_tmp=zero
 do iq1grad=1,nq1grad
   qdir=q1grad(2,iq1grad)
   do iatpert=1,natpert

     !Calculate the gradient of the potential generated by the first order atomic displacement density
     rhog1_tmp(:,:)=rhog1_atdis(iatpert,:,:)
     call hartredq(2,gmet,gsqcut,mpi_enreg,nfft,ngfft,qdir,rhog1_tmp,vqgradhart) 

     !To ckeck
     !call appdig(q2grad(3,iq2grad)+q1grad(3,iq1grad),"Gradient_Hartree_potential",fi1o)
     !call fftdatar_write_from_hdr("first_order_potential",fi1o,dtset%iomode,hdr_den,&
     ! & ngfft,cplex,nfft,nspden,vqgradhart,mpi_enreg)

     do iq2grad=1,nq2grad

       !Calculate the electrostatic energy term with the first order electric field density 
       if (timrev==1) then
         do ii=1,nfft
           jj=ii*2
           rhor1_tmp(jj-1,:)=rhor1_efield(iq2grad,ii,:)
         end do
       else if (timrev==0) then
         rhor1_tmp(:,:)=rhor1_efield(iq2grad,:,:)
       end if
       
       call dotprod_vn(2,rhor1_tmp,dotr,doti,nfft,nfftot,nspden,2,vqgradhart,ucvol)
       eqgradhart(re,iatpert,iq2grad,iq1grad)=dotr*half
       eqgradhart(im,iatpert,iq2grad,iq1grad)=doti*half

       qdrflg(pert_atdis(1,iatpert),pert_atdis(2,iatpert),q2grad(2,iq2grad),q1grad(2,iq1grad))=1

       blkflg(q2grad(2,iq2grad),q2grad(1,iq2grad),pert_atdis(2,iatpert),pert_atdis(1,iatpert),&
     &        q1grad(2,iq1grad),matom+8)=1           

     end do
   end do
 end do 

 ABI_DEALLOCATE(rhor1_tmp)
 ABI_DEALLOCATE(rhog1_tmp)
 ABI_DEALLOCATE(rhog1_atdis)
 ABI_DEALLOCATE(rhor1_efield)

!################# WAVE FUNCTION CONTRIBUTIONS  #######################################

!Determine the subset of symmetry operations (nsym1 operations)
!that leaves the perturbation invariant, and initialize corresponding arrays
!symaf1, symrl1, tnons1 (and pawang1%zarot, if PAW)..
!MR TODO: For the moment only the identiy symmetry is activated (nsym1=1) 
!         In a future I will try to activate perturbation dependent symmetries
!         with littlegroup_pert.F90. 
 nsym1 = 1
 ABI_ALLOCATE(indsy1,(4,nsym1,dtset%natom))
 ABI_ALLOCATE(symrc1,(3,3,nsym1))
 ABI_ALLOCATE(symaf1,(nsym1))
 ABI_ALLOCATE(symrl1,(3,3,nsym1))
 ABI_ALLOCATE(tnons1,(3,nsym1))
 symaf1(1:nsym1)= 1
 symrl1(:,:,nsym1)= dtset%symrel(:,:,1)
 tnons1(:,nsym1)= 0_dp

!Set up corresponding symmetry data
 ABI_ALLOCATE(irrzon1,(dtset%nfft**(1-1/nsym1),2,(nspden/dtset%nsppol)-3*(nspden/4)))
 ABI_ALLOCATE(phnons1,(2,dtset%nfft**(1-1/nsym1),(nspden/dtset%nsppol)-3*(nspden/4)))
 call setsym(indsy1,irrzon1,1,dtset%natom,dtset%nfft,dtset%ngfft,nspden,dtset%nsppol,&
&nsym1,phnons1,symaf1,symrc1,symrl1,tnons1,dtset%typat,xred)

 ABI_DEALLOCATE(indsy1)
 ABI_DEALLOCATE(symaf1)
 ABI_DEALLOCATE(symrl1)
 ABI_DEALLOCATE(tnons1)

!Determine the subset of k-points needed in the "reduced Brillouin zone",
!and initialize other quantities
 ABI_ALLOCATE(indkpt1_tmp,(nkpt))
 ABI_ALLOCATE(wtk_folded,(nkpt))
 ABI_ALLOCATE(bz2ibz_smap,(6, nkpt))
 indkpt1_tmp(:)=0 

 if (dtset%kptopt==2) then
   call symkpt(0,gmet,indkpt1_tmp,ab_out,dtset%kptns,nkpt,nkpt_rbz,&
  & nsym1,symrc1,timrev,dtset%wtk,wtk_folded,bz2ibz_smap,xmpi_comm_self)
 else if (dtset%kptopt==3) then
   call symkpt(0,gmet,indkpt1_tmp,ab_out,dtset%kptns,nkpt,nkpt_rbz,&
  & nsym1,symrc1,0,dtset%wtk,wtk_folded,bz2ibz_smap,xmpi_comm_self)
 else
   write(msg,"(1a)") 'kptopt must be 2 or 3 for the quadrupole calculation'
   MSG_BUG(msg)
 end if
 ABI_DEALLOCATE(bz2ibz_smap)

 ABI_ALLOCATE(doccde_rbz,(dtset%mband*nkpt_rbz*dtset%nsppol))
 ABI_ALLOCATE(indkpt1,(nkpt_rbz))
 ABI_ALLOCATE(istwfk_rbz,(nkpt_rbz))
 ABI_ALLOCATE(kpt_rbz,(3,nkpt_rbz))
 ABI_ALLOCATE(nband_rbz,(nkpt_rbz*dtset%nsppol))
 ABI_ALLOCATE(occ_rbz,(dtset%mband*nkpt_rbz*dtset%nsppol))
 ABI_ALLOCATE(wtk_rbz,(nkpt_rbz))
 indkpt1(:)=indkpt1_tmp(1:nkpt_rbz)
 do ikpt=1,nkpt_rbz
     istwfk_rbz(ikpt)=dtset%istwfk(indkpt1(ikpt))
     kpt_rbz(:,ikpt)=dtset%kptns(:,indkpt1(ikpt))
     wtk_rbz(ikpt)=wtk_folded(indkpt1(ikpt))
 end do
 ABI_DEALLOCATE(indkpt1_tmp)
 ABI_DEALLOCATE(wtk_folded)

!Transfer occ to occ_rbz 
!NOTE : this takes into account that indkpt1 is ordered
!MG: What about using occ(band,kpt,spin) ???
 bdtot_index=0;bdtot1_index=0
 do isppol=1,dtset%nsppol
   ikpt1=1
   do ikpt=1,nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
!    Must test against ikpt1/=nkpt_rbz+1, before evaluate indkpt1(ikpt1)
     if(ikpt1/=nkpt_rbz+1)then
       if(ikpt==indkpt1(ikpt1))then
         nband_rbz(ikpt1+(isppol-1)*nkpt_rbz)=nband_k
         occ_rbz(1+bdtot1_index:nband_k+bdtot1_index)=occ(1+bdtot_index:nband_k+bdtot_index)
         doccde_rbz(1+bdtot1_index:nband_k+bdtot1_index)=doccde(1+bdtot_index:nband_k+bdtot_index)
         ikpt1=ikpt1+1
         bdtot1_index=bdtot1_index+nband_k
       end if
     end if
     bdtot_index=bdtot_index+nband_k
   end do
 end do

!Compute maximum number of planewaves at k
call getmpw(ecut_eff,dtset%exchn2n3d,gmet,istwfk_rbz,kpt_rbz,mpi_enreg,mpw,nkpt_rbz)

!Allocate some k-dependent arrays at k
 ABI_ALLOCATE(kg,(3,mpw*nkpt_rbz))
 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(npwarr,(nkpt_rbz))
 ABI_ALLOCATE(npwtot,(nkpt_rbz))

!Determine distribution of k-points/bands over MPI processes
 if (allocated(mpi_enreg%my_kpttab)) then
   ABI_DEALLOCATE(mpi_enreg%my_kpttab)
 end if
 ABI_ALLOCATE(mpi_enreg%my_kpttab,(nkpt_rbz))
 if(xmpi_paral==1) then
   ABI_ALLOCATE(mpi_enreg%proc_distrb,(nkpt_rbz,dtset%mband,dtset%nsppol))
   call distrb2(dtset%mband,nband_rbz,nkpt_rbz,mpi_enreg%nproc_cell,dtset%nsppol,mpi_enreg)
 else
   mpi_enreg%my_kpttab(:)=(/(ii,ii=1,nkpt_rbz)/)
 end if
 my_nkpt_rbz=maxval(mpi_enreg%my_kpttab)
 call initmpi_band(mpi_enreg,nband_rbz,nkpt_rbz,dtset%nsppol)
 mkmem_rbz =my_nkpt_rbz 
 
!Set up the basis sphere of planewaves at k
 call kpgio(ecut_eff,dtset%exchn2n3d,gmet,istwfk_rbz,kg,&
& kpt_rbz,mkmem_rbz,nband_rbz,nkpt_rbz,'PERS',mpi_enreg,mpw,npwarr,npwtot,dtset%nsppol)
 ABI_DEALLOCATE(npwtot)

!Set up the spherical harmonics (Ylm) and 1st gradients at k
 useylmgr=1; option=1 ; nylmgr=3
 ABI_ALLOCATE(ylm,(mpw*mkmem_rbz,psps%mpsang*psps%mpsang*psps%useylm))
 ABI_ALLOCATE(ylmgr,(mpw*mkmem_rbz,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr))
 if (psps%useylm==1) then
   call initylmg(gprimd,kg,kpt_rbz,mkmem_rbz,mpi_enreg,psps%mpsang,mpw,nband_rbz,nkpt_rbz,&
&   npwarr,dtset%nsppol,option,rprimd,ylm,ylmgr)
 end if

!Initialize band structure datatype at k
 bantot_rbz=sum(nband_rbz(1:nkpt_rbz*dtset%nsppol))
 ABI_ALLOCATE(eigen0,(bantot_rbz))
 eigen0(:)=zero
 call ebands_init(bantot_rbz,bs_rbz,dtset%nelect,doccde_rbz,eigen0,istwfk_rbz,kpt_rbz,&
& nband_rbz,nkpt_rbz,npwarr,dtset%nsppol,dtset%nspinor,dtset%tphysel,dtset%tsmear,dtset%occopt,occ_rbz,wtk_rbz,&
& dtset%charge, dtset%kptopt, dtset%kptrlatt_orig, dtset%nshiftk_orig, dtset%shiftk_orig, &
& dtset%kptrlatt, dtset%nshiftk, dtset%shiftk)
 ABI_DEALLOCATE(eigen0)

 ABI_DEALLOCATE(doccde_rbz)

!Initialize header, update it with evolving variables
 gscase=0 ! A GS WF file is read
 call hdr_init(bs_rbz,codvsn,dtset,hdr0,pawtab,gscase,psps,wvl%descr,&
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

 call hdr0%update(bantot_rbz,etotal,fermie,&
& residm,rprimd,occ_rbz,pawrhoij,xred,dtset%amu_orig(:,1),&
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

!Clean band structure datatype (should use it more in the future !)
 call ebands_free(bs_rbz)

!!Initialize GS wavefunctions at k
 ireadwf0=1; formeig=0 ; ask_accurate=1 ; optorth=0
 mcg=mpw*dtset%nspinor*dtset%mband*mkmem_rbz*dtset%nsppol
 if (one*mpw*dtset%nspinor*dtset%mband*mkmem_rbz*dtset%nsppol > huge(1)) then
   write (msg,'(4a, 5(a,i0), 2a)')&
&   "Default integer is not wide enough to store the size of the GS wavefunction array (WF0, mcg).",ch10,&
&   "Action: increase the number of processors. Consider also OpenMP threads.",ch10,&
&   "nspinor: ",dtset%nspinor, "mpw: ",mpw, "mband: ",dtset%mband, "mkmem_rbz: ",&
&   mkmem_rbz, "nsppol: ",dtset%nsppol,ch10,&
&   'Note: Compiling with large int (int64) requires a full software stack (MPI/FFTW/BLAS/LAPACK...) compiled in int64 mode'
   MSG_ERROR(msg)
 end if
 ABI_STAT_ALLOCATE(cg,(2,mcg), ierr)
 ABI_CHECK(ierr==0, "out-of-memory in cg")

 ABI_ALLOCATE(eigen0,(dtset%mband*nkpt_rbz*dtset%nsppol))
 call inwffil(ask_accurate,cg,dtset,dtset%ecut,ecut_eff,eigen0,dtset%exchn2n3d,&
& formeig,hdr0,ireadwf0,istwfk_rbz,kg,&
& kpt_rbz,dtset%localrdwf,dtset%mband,mcg,&
& mkmem_rbz,mpi_enreg,mpw,nband_rbz,dtset%ngfft,nkpt_rbz,npwarr,&
& dtset%nsppol,dtset%nsym,occ_rbz,optorth,dtset%symafm,&
& dtset%symrel,dtset%tnons,dtfil%unkg,wffgs,wfftgs,&
& dtfil%unwffgs,dtfil%fnamewffk,wvl)
 ABI_DEALLOCATE(eigen0)

!Close wffgs%unwff, if it was ever opened (in inwffil)
 if (ireadwf0==1) then
   call WffClose(wffgs,ierr)
 end if

!==== Initialize most of the Hamiltonian ====
!1) Allocate all arrays and initialize quantities that do not depend on k and spin.
!2) Perform the setup needed for the non-local factors:
!3) Constant kleimann-Bylander energies are copied from psps to gs_hamkq.
 call init_hamiltonian(gs_hamkq,psps,pawtab,dtset%nspinor,nsppol,nspden,dtset%natom,&
& dtset%typat,xred,nfft,dtset%mgfft,ngfft,rprimd,dtset%nloalg,ph1d=ph1d,&
& use_gpu_cuda=dtset%use_gpu_cuda)


!==== Initialize response functions files and handlers ====
 !Atomic displacement files
 ABI_ALLOCATE(wfk_t_atdis,(natpert))

 do iatpert=1,natpert

   pertcase=pert_atdis(3,iatpert)
   call appdig(pertcase,dtfil%fnamewff1,fiwfatdis)

   !The value 20 is taken arbitrarily I would say
   forunit=20+pertcase

   !Check that atdis file exists and open it
   if (.not. file_exists(fiwfatdis)) then
     ! Trick needed to run Abinit test suite in netcdf mode.
     if (file_exists(nctk_ncify(fiwfatdis))) then
       write(msg,"(3a)")"- File: ",trim(fiwfatdis),&
       " does not exist but found netcdf file with similar name."
       call wrtout(std_out,msg,'COLL')
       fiwfatdis = nctk_ncify(fiwfatdis)
     end if
     if (.not. file_exists(fiwfatdis)) then
       MSG_ERROR('Missing file: '//TRIM(fiwfatdis))
     end if
   end if
   write(msg,'(a,a)')'-open atomic displacement wf1 file :',trim(fiwfatdis)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   call wfk_open_read(wfk_t_atdis(iatpert),fiwfatdis,formeig1,dtset%iomode,forunit,spaceworld)

 end do 


 !ddk files
 ABI_ALLOCATE(wfk_t_ddk,(nq1grad))
 do iq1grad=1,nq1grad

   pertcase=q1grad(3,iq1grad)
   call appdig(pertcase,dtfil%fnamewffddk,fiwfddk)

   !The value 20 is taken arbitrarily I would say
   forunit=20+pertcase

   !Check that ddk file exists and open it
   if (.not. file_exists(fiwfddk)) then
     ! Trick needed to run Abinit test suite in netcdf mode.
     if (file_exists(nctk_ncify(fiwfddk))) then
       write(msg,"(3a)")"- File: ",trim(fiwfddk),&
       " does not exist but found netcdf file with similar name."
       call wrtout(std_out,msg,'COLL')
       fiwfddk = nctk_ncify(fiwfddk)
     end if
     if (.not. file_exists(fiwfddk)) then
       MSG_ERROR('Missing file: '//TRIM(fiwfddk))
     end if
   end if
   write(msg, '(a,a)') '-open ddk wf1 file :',trim(fiwfddk)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   call wfk_open_read(wfk_t_ddk(iq1grad),fiwfddk,formeig1,dtset%iomode,forunit,spaceworld)

 end do 

 !Electric field files
 ABI_ALLOCATE(wfk_t_efield,(nq2grad))
 do iq2grad=1,nq2grad

   pertcase=q2grad(3,iq2grad)
   call appdig(pertcase,dtfil%fnamewff1,fiwfefield)

   !The value 20 is taken arbitrarily I would say
   forunit=20+pertcase

   !Check that efield file exists and open it
   if (.not. file_exists(fiwfefield)) then
     ! Trick needed to run Abinit test suite in netcdf mode.
     if (file_exists(nctk_ncify(fiwfefield))) then
       write(msg,"(3a)")"- File: ",trim(fiwfefield),&
       " does not exist but found netcdf file with similar name."
       call wrtout(std_out,msg,'COLL')
       fiwfefield = nctk_ncify(fiwfefield)
     end if
     if (.not. file_exists(fiwfefield)) then
       MSG_ERROR('Missing file: '//TRIM(fiwfefield))
     end if
   end if
   write(msg, '(a,a)') '-open electric field wf1 file :',trim(fiwfefield)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   call wfk_open_read(wfk_t_efield(iq2grad),fiwfefield,formeig1,dtset%iomode,forunit,spaceworld)

 end do

 !d2_dkdk
 ABI_ALLOCATE(wfk_t_dkdk,(nq1q2grad))
 do iq1q2grad=1,nq1q2grad

   pertcase=q1q2grad(4,iq1q2grad)
   call appdig(pertcase,dtfil%fnamewffdkdk,fiwfdkdk)

   !The value 20 is taken arbitrarily I would say
   forunit=20+pertcase

   !Check that d2_ddk file exists and open it
   if (.not. file_exists(fiwfdkdk)) then
     ! Trick needed to run Abinit test suite in netcdf mode.
     if (file_exists(nctk_ncify(fiwfdkdk))) then
       write(msg,"(3a)")"- File: ",trim(fiwfdkdk),&
       " does not exist but found netcdf file with similar name."
       call wrtout(std_out,msg,'COLL')
       fiwfdkdk = nctk_ncify(fiwfdkdk)
     end if
     if (.not. file_exists(fiwfdkdk)) then
       MSG_ERROR('Missing file: '//TRIM(fiwfdkdk))
     end if
   end if
   if (iq1q2grad <= 6) then
     write(msg, '(a,a)') '-open d2_dkdk wf2 file :',trim(fiwfdkdk)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
     call wfk_open_read(wfk_t_dkdk(iq1q2grad),fiwfdkdk,formeig1,dtset%iomode,forunit,spaceworld)
   else
     wfk_t_dkdk(iq1q2grad)=wfk_t_dkdk(iq1q2grad-3)
   end if

 end do

!Allocate the quadrupole tensor part depending on the wave functions
 ABI_ALLOCATE(qdrpwf,(2,natpert,nq2grad,nq1grad))
 ABI_ALLOCATE(qdrpwf_k,(2,natpert,nq2grad,nq1grad))
 ABI_ALLOCATE(qdrpwf_t1,(2,natpert,nq2grad,nq1grad))
 ABI_ALLOCATE(qdrpwf_t1_k,(2,natpert,nq2grad,nq1grad))
 ABI_ALLOCATE(qdrpwf_t2,(2,natpert,nq2grad,nq1grad))
 ABI_ALLOCATE(qdrpwf_t2_k,(2,natpert,nq2grad,nq1grad))
 ABI_ALLOCATE(qdrpwf_t3,(2,natpert,nq2grad,nq1grad))
 ABI_ALLOCATE(qdrpwf_t3_k,(2,natpert,nq2grad,nq1grad))
 ABI_ALLOCATE(qdrpwf_t4,(2,natpert,nq2grad,nq1grad))
 ABI_ALLOCATE(qdrpwf_t4_k,(2,natpert,nq2grad,nq1grad))
 ABI_ALLOCATE(qdrpwf_t5,(2,natpert,nq2grad,nq1grad))
 ABI_ALLOCATE(qdrpwf_t5_k,(2,natpert,nq2grad,nq1grad))
 qdrpwf=zero
 qdrpwf_t1=zero
 qdrpwf_t2=zero
 qdrpwf_t3=zero
 qdrpwf_t4=zero
 qdrpwf_t5=zero

!LOOP OVER SPINS
 bdtot_index=0
 icg=0
 do isppol=1,nsppol
   ikg=0

!  Continue to initialize the Hamiltonian
   call gs_hamkq%load_spin(isppol,with_nonlocal=.true.)

!  BIG FAT k POINT LOOP
   do ikpt=1,nkpt_rbz

     nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
     istwf_k=istwfk_rbz(ikpt)
     npw_k=npwarr(ikpt)
     npw1_k=npw_k

     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
       bdtot_index=bdtot_index+nband_k

       cycle ! Skip the rest of the k-point loop
     end if

     ABI_ALLOCATE(occ_k,(nband_k))
     ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
     ABI_ALLOCATE(ylmgr_k,(npw_k,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr))
     occ_k(:)=occ_rbz(1+bdtot_index:nband_k+bdtot_index)
     kpoint(:)=kpt_rbz(:,ikpt)
     wtk_k=wtk_rbz(ikpt)

!    Get plane-wave vectors and related data at k
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     if (psps%useylm==1) then
       do ilm=1,psps%mpsang*psps%mpsang
         ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       end do
       if (useylmgr==1) then
         do ilm=1,psps%mpsang*psps%mpsang
           do ii=1,nylmgr
             ylmgr_k(1:npw_k,ii,ilm)=ylmgr(1+ikg:npw_k+ikg,ii,ilm)
           end do
         end do
       end if
     end if

     call dfpt_qdrpwf(atindx,cg,cplex,dtset,gs_hamkq,gsqcut, &
     &  icg,ikpt,indkpt1,isppol,istwf_k, &
     &  kg_k,kpoint,mkmem_rbz, &
     &  mpi_enreg,mpw,natpert,nattyp,nband_k,nfft,ngfft,nkpt_rbz, &
     &  npw_k,nq1grad, &
     &  nq2grad,nq1q2grad,nspden,nsppol,nylmgr,occ_k, &
     &  pert_atdis,ph1d,psps,qdrpwf_k,qdrpwf_t1_k,qdrpwf_t2_k,qdrpwf_t3_k,    &
     &  qdrpwf_t4_k,qdrpwf_t5_k,q1grad,q2grad,q1q2grad,rmet,ucvol,useylmgr, &
     &  vhxc1_atdis,vhxc1_efield,wfk_t_atdis,wfk_t_efield,wfk_t_ddk, &
     &  wfk_t_dkdk,wtk_k,xred,ylm_k,ylmgr_k)


!    Add the contribution from each k-point
     qdrpwf=qdrpwf + qdrpwf_k
     qdrpwf_t1=qdrpwf_t1 + qdrpwf_t1_k
     qdrpwf_t2=qdrpwf_t2 + qdrpwf_t2_k
     qdrpwf_t3=qdrpwf_t3 + qdrpwf_t3_k
     qdrpwf_t4=qdrpwf_t4 + qdrpwf_t4_k
     qdrpwf_t5=qdrpwf_t5 + qdrpwf_t5_k

!    Keep track of total number of bands
     bdtot_index=bdtot_index+nband_k

!    Shift arrays memory
     if (mkmem_rbz/=0) then
       icg=icg+npw_k*dtset%nspinor*nband_k
       ikg=ikg+npw_k
     end if

     ABI_DEALLOCATE(occ_k)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylmgr_k)

   end do
!  END BIG FAT k POINT LOOP
 end do
!END LOOP OVER SPINS


!Close response function files
 do iatpert=1,natpert
   call wfk_t_atdis(iatpert)%close()
 end do
 do iq1grad=1,nq1grad
   call wfk_t_ddk(iq1grad)%close()
 end do
 do iq2grad=1,nq2grad
   call wfk_t_efield(iq2grad)%close()
 end do
 do iq1q2grad=1,nq1q2grad
   if (iq1q2grad <= 6) call wfk_t_dkdk(iq1q2grad)%close()
 end do

!=== MPI communications ==================
 if (xmpi_paral==1) then

   call xmpi_sum(qdrpwf,spaceworld,ierr)
   call xmpi_sum(qdrpwf_t1,spaceworld,ierr)
   call xmpi_sum(qdrpwf_t2,spaceworld,ierr)
   call xmpi_sum(qdrpwf_t3,spaceworld,ierr)
   call xmpi_sum(qdrpwf_t4,spaceworld,ierr)
   call xmpi_sum(qdrpwf_t5,spaceworld,ierr)

 end if

!Anounce finalization of quadrupole tensor calculation
 write(msg, '(a,a,a)' ) ch10, &
' Quadrupole tensor calculation completed ',ch10
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

!Gather the different terms in the quadrupole tensor and print them out
 if (me==0) then
 call dfpt_qdrpout(d3etot,eqgradhart,gprimd,dtset%kptopt,matom,mpert,natpert,& 
    & nq1grad,nq2grad,pert_atdis,dtset%prtvol,q1grad,q2grad,qdrflg,qdrpwf,qdrpwf_t1,qdrpwf_t2, &
    & qdrpwf_t3,qdrpwf_t4,qdrpwf_t5,rprimd,ucvol)
 end if

 !Deallocations
 ABI_DEALLOCATE(eqgradhart)
 ABI_DEALLOCATE(indkpt1)
 ABI_DEALLOCATE(istwfk_rbz)
 ABI_DEALLOCATE(kg)
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(kpt_rbz)
! ABI_DEALLOCATE(mpi_enreg%my_kpttab)
! ABI_DEALLOCATE(mpi_enreg%proc_distrb)
 ABI_DEALLOCATE(nband_rbz)
 ABI_DEALLOCATE(npwarr)
 ABI_DEALLOCATE(occ_rbz)
 ABI_DEALLOCATE(ph1d)
 ABI_DEALLOCATE(pert_atdis)
 ABI_DEALLOCATE(qdrflg)
 ABI_DEALLOCATE(qdrpwf)
 ABI_DEALLOCATE(qdrpwf_k)
 ABI_DEALLOCATE(qdrpwf_t1)
 ABI_DEALLOCATE(qdrpwf_t2)
 ABI_DEALLOCATE(qdrpwf_t3)
 ABI_DEALLOCATE(qdrpwf_t4)
 ABI_DEALLOCATE(qdrpwf_t5)
 ABI_DEALLOCATE(qdrpwf_t1_k)
 ABI_DEALLOCATE(qdrpwf_t2_k)
 ABI_DEALLOCATE(qdrpwf_t3_k)
 ABI_DEALLOCATE(qdrpwf_t4_k)
 ABI_DEALLOCATE(qdrpwf_t5_k)
 ABI_DEALLOCATE(q1grad)
 ABI_DEALLOCATE(q1q2grad)
 ABI_DEALLOCATE(q2grad)
 ABI_DEALLOCATE(vhxc1_atdis)
 ABI_DEALLOCATE(vhxc1_efield)
 ABI_DEALLOCATE(vqgradhart)
 ABI_DEALLOCATE(wfk_t_atdis)
 ABI_DEALLOCATE(wfk_t_ddk)
 ABI_DEALLOCATE(wfk_t_dkdk)
 ABI_DEALLOCATE(wfk_t_efield)
 ABI_DEALLOCATE(wtk_rbz)
 ABI_DEALLOCATE(ylm)
 ABI_DEALLOCATE(ylmgr)
 if(xmpi_paral==1) then
   ABI_DEALLOCATE(mpi_enreg%proc_distrb)
 end if

 DBG_EXIT("COLL")

end subroutine dfpt_qdrpole
!!***

!!****f* ABINIT/dfpt_qdrpout
!! NAME
!!  dfpt_qdrpout
!!
!! FUNCTION
!!  This subroutine gathers the different terms entering the quadrupole tensor,
!!  perfofms the transformation from reduced to cartesian coordinates and 
!!  writes out the quadrupole tensor in external files. It also writes the firts
!!  q-gradient of the polarization response to an atomic displacement. 
!!  
!! COPYRIGHT
!!  Copyright (C) 2018 ABINIT group (MR,MS)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  eqgradhart(2,natpert,nq2grad,nq1grad)=electrostatic contribution from the 
!!                                             q-gradient of the Hartree potential
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  kptopt=2 time reversal symmetry is enforced, 3 trs is not enforced (for debugging purposes)
!!  matom=maximum number of atoms to displace
!!  natpert=number of atomic displacement perturbations
!!  nq1grad=number of q1 (q_{\gamma}) gradients
!!  nq2grad=number of q2 (q_{\delta}) gradients
!!  pert_atdis(3,natpert)=array with the info for the atomic displacement perturbations
!!  prtvol=volume of information to be printed. 1-> The different contributions to the quadrupole are printed.
!!  q1grad(3,nq1grad)=array with the info for the q1 (q_{\gamma}) gradients
!!  q2grad(3,nq2grad)=array with the info for the q2 (q_{\gamma}) gradients
!!  qdrflg(matom,3,3,3)=array that indicates which quadrupole tensor elements have been calculated
!!  qdrpwf_k(2,natpert,nq2grad,nq1grad)= wave function dependent part of the quadrupole tensor
!!                                       for the k-point kpt
!!  qdrpwf_t1_k(2,natpert,nq2grad,nq1grad)= t1 term (see notes) of qdrpwf_k
!!  qdrpwf_t2_k(2,natpert,nq2grad,nq1grad)= t2 term (see notes) of qdrpwf_k
!!  qdrpwf_t3_k(2,natpert,nq2grad,nq1grad)= t3 term (see notes) of qdrpwf_k
!!  qdrpwf_t4_k(2,natpert,nq2grad,nq1grad)= t4 term (see notes) of qdrpwf_k
!!  qdrpwf_t5_k(2,natpert,nq2grad,nq1grad)= t5 term (see notes) of qdrpwf_k
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  ucvol=unit cell volume in bohr**3.
!!  
!! OUTPUT
!!  d3etot(2,3,mpert,6,mpert,3,mpert)= matrix of the 3DTE
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!!  dfpt_qdrpole
!!
!! CHILDREN
!!
!!  cart39
!!
!! SOURCE

subroutine dfpt_qdrpout(d3etot,eqgradhart,gprimd,kptopt,matom,mpert,natpert, & 
         & nq1grad,nq2grad,pert_atdis,prtvol,q1grad,q2grad,qdrflg,qdrpwf,qdrpwf_t1,qdrpwf_t2, & 
         & qdrpwf_t3,qdrpwf_t4,qdrpwf_t5,rprimd,ucvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_qdrpout'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: kptopt,matom,mpert,natpert,nq1grad,nq2grad,prtvol
 real(dp),intent(inout) :: d3etot(2,3,mpert,3,mpert,3,mpert)
 real(dp),intent(in) :: ucvol

!arrays
 integer,intent(in) :: pert_atdis(3,natpert)
 integer,intent(in) :: qdrflg(matom,3,3,3)
 integer,intent(in) :: q1grad(3,nq1grad),q2grad(3,nq2grad)
 real(dp),intent(in) :: eqgradhart(2,natpert,nq2grad,nq1grad)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: qdrpwf(2,natpert,nq2grad,nq1grad)
 real(dp),intent(in) :: qdrpwf_t1(2,natpert,nq2grad,nq1grad)
 real(dp),intent(in) :: qdrpwf_t2(2,natpert,nq2grad,nq1grad)
 real(dp),intent(in) :: qdrpwf_t3(2,natpert,nq2grad,nq1grad)
 real(dp),intent(in) :: qdrpwf_t4(2,natpert,nq2grad,nq1grad)
 real(dp),intent(in) :: qdrpwf_t5(2,natpert,nq2grad,nq1grad)
 real(dp),intent(in) :: rprimd(3,3)
 
!Local variables-------------------------------
!scalars
 integer :: alpha,beta,gamma
 integer :: iatpert,iatdir,iatom,ii,iiq1grad,iiq2grad
 integer :: iq1dir,iq1grad,iq2dir,iq2grad 
 integer :: iq1pert,iq2pert
 integer, parameter :: re=1,im=2
 real(dp) :: piezore,piezoim,tmpre, tmpim
 character(len=500) :: msg
!arrays
 integer :: flg1(3),flg2(3)
 integer,allocatable :: cartflg(:,:,:,:)
 real(dp) :: vec1(3),vec2(3)
 real(dp),allocatable :: qdrptens_cart(:,:,:,:,:),qdrptens_red(:,:,:,:,:)
 real(dp),allocatable :: dqpol_cart(:,:,:,:,:),dqpol_red(:,:,:,:,:)
 
! *************************************************************************

 DBG_ENTER("COLL")

!Open output files
 if (prtvol==1) then
   open(unit=71,file='qdrpl_wf_t1.out',status='unknown',form='formatted',action='write')
   open(unit=72,file='qdrpl_wf_t2.out',status='unknown',form='formatted',action='write')
   open(unit=73,file='qdrpl_wf_t3.out',status='unknown',form='formatted',action='write')
   open(unit=74,file='qdrpl_wf_t4.out',status='unknown',form='formatted',action='write')
   open(unit=75,file='qdrpl_wf_t5.out',status='unknown',form='formatted',action='write')
   open(unit=76,file='qdrpl_elecstic.out',status='unknown',form='formatted',action='write')
 end if

!Gather the different terms in the tensors and print the result
 ABI_ALLOCATE(qdrptens_red,(2,matom,3,3,3))
 ABI_ALLOCATE(dqpol_red,(2,matom,3,3,3))

 if (kptopt==3) then

   !Write real and 'true' imaginary parts of quadrupole tensor and the 
   !q-gradient of the polarization response
   iq1pert=matom+8
   iq2pert=matom+2
   do iq1grad=1,nq1grad
     iq1dir=q1grad(2,iq1grad)
     do iq2grad=1,nq2grad
       iq2dir=q2grad(2,iq2grad)
       do iatpert=1,natpert
         iatom=pert_atdis(1,iatpert)
         iatdir=pert_atdis(2,iatpert)

         if (qdrflg(iatom,iatdir,iq2dir,iq1dir)==1) then

           !Calculate and save the third order energy derivative
           tmpre=eqgradhart(re,iatpert,iq2grad,iq1grad)+qdrpwf(re,iatpert,iq2grad,iq1grad)
           tmpim=eqgradhart(im,iatpert,iq2grad,iq1grad)+qdrpwf(im,iatpert,iq2grad,iq1grad)
           d3etot(1,iq2dir,iq2pert,iatdir,iatom,iq1dir,iq1pert)=tmpre
           d3etot(2,iq2dir,iq2pert,iatdir,iatom,iq1dir,iq1pert)=tmpim

           !Calculate and write the q-gradient of the polarization response 
           !(without the inverse volume factor)
           dqpol_red(1,iatom,iatdir,iq2dir,iq1dir)=-two*tmpim
           dqpol_red(2,iatom,iatdir,iq2dir,iq1dir)=two*tmpre
 
           if (qdrflg(iatom,iatdir,iq2dir,iq1dir)==1 .and. qdrflg(iatom,iatdir,iq1dir,iq2dir)==1 ) then

             !Avoid double counting when summing up unsymmetrized contributions
             do iiq1grad=1,3
               if (q2grad(2,iq2grad)==q1grad(2,iiq1grad)) exit
             end do

             do iiq2grad=1,3
               if (q1grad(2,iq1grad)==q2grad(2,iiq2grad)) exit
             end do

             qdrptens_red(re,iatom,iatdir,iq2dir,iq1dir)=-two* &
           & ( eqgradhart(re,iatpert,iq2grad,iq1grad)+eqgradhart(re,iatpert,iiq2grad,iiq1grad) &
           & + qdrpwf(re,iatpert,iq2grad,iq1grad)+qdrpwf(re,iatpert,iiq2grad,iiq1grad) ) 
             qdrptens_red(im,iatom,iatdir,iq2dir,iq1dir)=-two* &
           & ( eqgradhart(im,iatpert,iq2grad,iq1grad)+eqgradhart(im,iatpert,iiq2grad,iiq1grad) &
           & + qdrpwf(im,iatpert,iq2grad,iq1grad)+qdrpwf(im,iatpert,iiq2grad,iiq1grad) ) 

             !Divide by the imaginary unit to get the proper observable 
             !(See M.Royo and M.Stengel paper)
             tmpre=qdrptens_red(re,iatom,iatdir,iq2dir,iq1dir)
             tmpim=qdrptens_red(im,iatom,iatdir,iq2dir,iq1dir)
             qdrptens_red(re,iatom,iatdir,iq2dir,iq1dir)=tmpim
             qdrptens_red(im,iatom,iatdir,iq2dir,iq1dir)=-tmpre

             if (prtvol==1) then
               !Write individual contributions 
               write(71,'(4(i2,2x),2(f18.10,2x))') iq2dir,iatom,iatdir,iq1dir,                  &
             & qdrpwf_t1(im,iatpert,iq2grad,iq1grad)+qdrpwf_t1(im,iatpert,iiq2grad,iiq1grad),   &
             & -(qdrpwf_t1(re,iatpert,iq2grad,iq1grad)+qdrpwf_t1(re,iatpert,iiq2grad,iiq1grad))

               write(72,'(4(i2,2x),2(f18.10,2x))') iq2dir,iatom,iatdir,iq1dir,                  &
             & qdrpwf_t2(im,iatpert,iq2grad,iq1grad)+qdrpwf_t2(im,iatpert,iiq2grad,iiq1grad),   &
             & -(qdrpwf_t2(re,iatpert,iq2grad,iq1grad)+qdrpwf_t2(re,iatpert,iiq2grad,iiq1grad))

               write(73,'(4(i2,2x),2(f18.10,2x))') iq2dir,iatom,iatdir,iq1dir,                  &
             & qdrpwf_t3(im,iatpert,iq2grad,iq1grad)+qdrpwf_t3(im,iatpert,iiq2grad,iiq1grad),   &
             & -(qdrpwf_t3(re,iatpert,iq2grad,iq1grad)+qdrpwf_t3(re,iatpert,iiq2grad,iiq1grad))

               write(74,'(4(i2,2x),2(f18.10,2x))') iq2dir,iatom,iatdir,iq1dir,                  &
             & qdrpwf_t4(im,iatpert,iq2grad,iq1grad)+qdrpwf_t4(im,iatpert,iiq2grad,iiq1grad),   &
             & -(qdrpwf_t4(re,iatpert,iq2grad,iq1grad)+qdrpwf_t4(re,iatpert,iiq2grad,iiq1grad))

               write(75,'(4(i2,2x),2(f18.10,2x))') iq2dir,iatom,iatdir,iq1dir,                  &
             & qdrpwf_t5(im,iatpert,iq2grad,iq1grad)+qdrpwf_t5(im,iatpert,iiq2grad,iiq1grad),   &
             & -(qdrpwf_t5(re,iatpert,iq2grad,iq1grad)+qdrpwf_t5(re,iatpert,iiq2grad,iiq1grad))

               write(76,'(4(i2,2x),2(f18.10,2x))') iq2dir,iatom,iatdir,iq1dir,                  &
             & eqgradhart(im,iatpert,iq2grad,iq1grad)+eqgradhart(im,iatpert,iiq2grad,iiq1grad), &
             & -(eqgradhart(re,iatpert,iq2grad,iq1grad)+eqgradhart(re,iatpert,iiq2grad,iiq1grad))
             end if

           end if
         end if
       end do
     end do
   end do

 else if (kptopt==2) then

   !Write real and zero imaginary parts of quadrupole tensor and the 
   !q-gradient of the polarization response
   iq1pert=matom+8
   iq2pert=matom+2
   do iq1grad=1,nq1grad
     iq1dir=q1grad(2,iq1grad)
     do iq2grad=1,nq2grad
       iq2dir=q2grad(2,iq2grad)
       do iatpert=1,natpert
         iatom=pert_atdis(1,iatpert)
         iatdir=pert_atdis(2,iatpert)

         if (qdrflg(iatom,iatdir,iq2dir,iq1dir)==1) then

           !Calculate and save the third order energy derivative
           tmpim=eqgradhart(im,iatpert,iq2grad,iq1grad)+qdrpwf(im,iatpert,iq2grad,iq1grad)
           d3etot(1,iq2dir,iq2pert,iatdir,iatom,iq1dir,iq1pert)=zero
           d3etot(2,iq2dir,iq2pert,iatdir,iatom,iq1dir,iq1pert)=tmpim

           !Calculate and write the q-gradient of the polarization response
           !(without the inverse volume factor)
           dqpol_red(1,iatom,iatdir,iq2dir,iq1dir)=-two*tmpim
           dqpol_red(2,iatom,iatdir,iq2dir,iq1dir)=0.0_dp

           if (qdrflg(iatom,iatdir,iq2dir,iq1dir)==1 .and. qdrflg(iatom,iatdir,iq1dir,iq2dir)==1 ) then

             do iiq1grad=1,3
               if (q2grad(2,iq2grad)==q1grad(2,iiq1grad)) exit
             end do

             do iiq2grad=1,3
               if (q1grad(2,iq1grad)==q2grad(2,iiq2grad)) exit
             end do

             qdrptens_red(re,iatom,iatdir,iq2dir,iq1dir)=-two* &
           & ( eqgradhart(re,iatpert,iq2grad,iq1grad)+eqgradhart(re,iatpert,iiq2grad,iiq1grad) &
           & + qdrpwf(re,iatpert,iq2grad,iq1grad)+qdrpwf(re,iatpert,iiq2grad,iiq1grad) )
             qdrptens_red(im,iatom,iatdir,iq2dir,iq1dir)=-two* &
           & ( eqgradhart(im,iatpert,iq2grad,iq1grad)+eqgradhart(im,iatpert,iiq2grad,iiq1grad) &
           & + qdrpwf(im,iatpert,iq2grad,iq1grad)+qdrpwf(im,iatpert,iiq2grad,iiq1grad) )

             !Divide by the imaginary unit to get the proper observable 
             !(See M.Royo and M.Stengel paper)
             tmpim=qdrptens_red(im,iatom,iatdir,iq2dir,iq1dir)
             qdrptens_red(re,iatom,iatdir,iq2dir,iq1dir)=tmpim
             qdrptens_red(im,iatom,iatdir,iq2dir,iq1dir)=0.0_dp

             if (prtvol==1) then
               !Write individual contributions 
               write(71,'(4(i2,2x),2(f18.10,2x))') iq2dir,iatom,iatdir,iq1dir,                  &
             & qdrpwf_t1(im,iatpert,iq2grad,iq1grad)+qdrpwf_t1(im,iatpert,iiq2grad,iiq1grad),   &
             & 0.0_dp

               write(72,'(4(i2,2x),2(f18.10,2x))') iq2dir,iatom,iatdir,iq1dir,                  & 
             & qdrpwf_t2(im,iatpert,iq2grad,iq1grad)+qdrpwf_t2(im,iatpert,iiq2grad,iiq1grad),   &
             & 0.0_dp

               write(73,'(4(i2,2x),2(f18.10,2x))') iq2dir,iatom,iatdir,iq1dir,                  &
             & qdrpwf_t3(im,iatpert,iq2grad,iq1grad)+qdrpwf_t3(im,iatpert,iiq2grad,iiq1grad),   &
             & 0.0_dp

               write(74,'(4(i2,2x),2(f18.10,2x))') iq2dir,iatom,iatdir,iq1dir,                  &
             & qdrpwf_t4(im,iatpert,iq2grad,iq1grad)+qdrpwf_t4(im,iatpert,iiq2grad,iiq1grad),   &
             & 0.0_dp

               write(75,'(4(i2,2x),2(f18.10,2x))') iq2dir,iatom,iatdir,iq1dir,                  &
             & qdrpwf_t5(im,iatpert,iq2grad,iq1grad)+qdrpwf_t5(im,iatpert,iiq2grad,iiq1grad),   &
             & 0.0_dp

               write(76,'(4(i2,2x),2(f18.10,2x))') iq2dir,iatom,iatdir,iq1dir,                  &
             & eqgradhart(im,iatpert,iq2grad,iq1grad)+eqgradhart(im,iatpert,iiq2grad,iiq1grad), &
             & 0.0_dp
             end if

           end if
         end if
       end do
     end do
   end do

 else
   write(msg,"(1a)") 'kptopt must be 2 or 3 for the quadrupole calculation'
   MSG_BUG(msg)
 end if

 if (prtvol==1) then
   close(71)
   close(72)
   close(73)
   close(74)
   close(75)
   close(76)
 end if

!Transformation to cartesian coordinates of the quadrupole tensor
 ABI_ALLOCATE(qdrptens_cart,(2,matom,3,3,3))
 ABI_ALLOCATE(dqpol_cart,(2,matom,3,3,3))
 ABI_ALLOCATE(cartflg,(matom,3,3,3))
 qdrptens_cart(:,:,:,:,:)=qdrptens_red(:,:,:,:,:)
 dqpol_cart(:,:,:,:,:)=dqpol_red(:,:,:,:,:)
 cartflg=0

 ABI_DEALLOCATE(qdrptens_red)
 ABI_DEALLOCATE(dqpol_red)

!1st transform coordenates of the atomic displacement derivative
 do iq1dir=1,3
   do iq2dir=1,3
     do ii=1,2
       do iatom=1,matom

         do iatdir=1,3
           vec1(iatdir)=qdrptens_cart(ii,iatom,iatdir,iq2dir,iq1dir)
           flg1(iatdir)=qdrflg(iatom,iatdir,iq2dir,iq1dir)
         end do 
         call cart39(flg1,flg2,gprimd,iatom,matom,rprimd,vec1,vec2)
         do iatdir=1,3
           qdrptens_cart(ii,iatom,iatdir,iq2dir,iq1dir)=vec2(iatdir)
           cartflg(iatom,iatdir,iq2dir,iq1dir)=flg2(iatdir)
         end do

         do iatdir=1,3
           vec1(iatdir)=dqpol_cart(ii,iatom,iatdir,iq2dir,iq1dir)
           flg1(iatdir)=qdrflg(iatom,iatdir,iq2dir,iq1dir)
         end do 
         call cart39(flg1,flg2,gprimd,iatom,matom,rprimd,vec1,vec2)
         do iatdir=1,3
           dqpol_cart(ii,iatom,iatdir,iq2dir,iq1dir)=vec2(iatdir)
         end do

       end do
     end do
   end do
 end do

!2nd transform coordinates of the electric field derivative
 do iq1dir=1,3
   do iatdir=1,3
     do iatom=1,matom
       do ii=1,2
         do iq2dir=1,3
           vec1(iq2dir)=qdrptens_cart(ii,iatom,iatdir,iq2dir,iq1dir)
           flg1(iq2dir)=qdrflg(iatom,iatdir,iq2dir,iq1dir)
         end do 
         call cart39(flg1,flg2,gprimd,matom+2,matom,rprimd,vec1,vec2)
         do iq2dir=1,3
           qdrptens_cart(ii,iatom,iatdir,iq2dir,iq1dir)=vec2(iq2dir)
           cartflg(iatom,iatdir,iq2dir,iq1dir)=flg2(iq2dir)
         end do

         do iq2dir=1,3
           vec1(iq2dir)=dqpol_cart(ii,iatom,iatdir,iq2dir,iq1dir)
           flg1(iq2dir)=qdrflg(iatom,iatdir,iq2dir,iq1dir)
         end do 
         call cart39(flg1,flg2,gprimd,matom+2,matom,rprimd,vec1,vec2)
         do iq2dir=1,3
           dqpol_cart(ii,iatom,iatdir,iq2dir,iq1dir)=vec2(iq2dir)
         end do

       end do
     end do
   end do
 end do

!3rd transform coordinates of the q-gradient (treat it as electric field)
 do iq2dir=1,3
   do iatdir=1,3
     do iatom=1,matom
       do ii=1,2
         do iq1dir=1,3
           vec1(iq1dir)=qdrptens_cart(ii,iatom,iatdir,iq2dir,iq1dir)
           flg1(iq1dir)=qdrflg(iatom,iatdir,iq2dir,iq1dir)
         end do
         call cart39(flg1,flg2,gprimd,matom+2,matom,rprimd,vec1,vec2)
         do iq1dir=1,3
           qdrptens_cart(ii,iatom,iatdir,iq2dir,iq1dir)=vec2(iq1dir)
           cartflg(iatom,iatdir,iq2dir,iq1dir)=flg2(iq1dir)
         end do

         do iq1dir=1,3
           vec1(iq1dir)=dqpol_cart(ii,iatom,iatdir,iq2dir,iq1dir)
           flg1(iq1dir)=qdrflg(iatom,iatdir,iq2dir,iq1dir)
         end do
         call cart39(flg1,flg2,gprimd,matom+2,matom,rprimd,vec1,vec2)
         do iq1dir=1,3
           dqpol_cart(ii,iatom,iatdir,iq2dir,iq1dir)=vec2(iq1dir)
         end do

       end do
     end do
   end do
 end do


!Write the Quadrupole tensor in cartesian coordinates
 write(ab_out,'(a)')' '
 write(ab_out,'(a)')' Quadrupole tensor, in cartesian coordinates,'
 write(ab_out,'(a)')' efidir   atom   atddir   qgrdir          real part        imaginary part'
 do iq1dir=1,3
   do iq2dir=1,3
     do iatdir=1,3
       do iatom=1,matom
        
         if (cartflg(iatom,iatdir,iq2dir,iq1dir)==1) then

           write(ab_out,'(4(i5,3x),2(1x,f20.10))') iq2dir,iatom,iatdir,iq1dir,                   &
         & qdrptens_cart(re,iatom,iatdir,iq2dir,iq1dir),qdrptens_cart(im,iatom,iatdir,iq2dir,iq1dir)

         end if

       end do
     end do
   end do
   write(ab_out,'(a)')' '
 end do

!Write the q-gradient of the Polarization response
 write(ab_out,*)' q-gradient of the polarization response '
 write(ab_out,*)' to an atomic displacementatom, in cartesian coordinates,'
 write(ab_out,*)' (1/ucvol factor not included),'
 write(ab_out,*)' efidir   atom   atddir   qgrdir          real part        imaginary part'
 do iq1dir=1,3
   do iq2dir=1,3
     do iatdir=1,3
       do iatom=1,matom
        
         if (cartflg(iatom,iatdir,iq2dir,iq1dir)==1) then

           write(ab_out,'(4(i5,3x),2(1x,f20.10))') iq2dir,iatom,iatdir,iq1dir,                  &
         & dqpol_cart(re,iatom,iatdir,iq2dir,iq1dir),dqpol_cart(im,iatom,iatdir,iq2dir,iq1dir)

         end if

       end do
     end do
   end do
   write(ab_out,*)' '
 end do

!Write the electronic (frozen-ion) contribution to the piezoelectric tensor
!(R.M. Martin, PRB 5, 1607 (1972))
 write(ab_out,'(a)')' '
 write(ab_out,'(a)')' Electronic (frozen-ion) contribution to the piezoelectric tensor,'
 write(ab_out,'(a)')' in cartesian coordinates, (from quadrupole calculation)'
 write(ab_out,'(a)')' atddir   qgrdir   efidir        real part           imaginary part'
 do iq1dir=1,3
   gamma=iq1dir
   do iq2dir=1,3
     alpha=iq2dir
     do iatdir=1,3
       beta=iatdir

       piezore=zero
       piezoim=zero
       do iatom=1,matom

         if (cartflg(iatom,iatdir,iq2dir,iq1dir)==1) then

           piezore=piezore+( qdrptens_cart(re,iatom,beta,alpha,gamma)-qdrptens_cart(re,iatom,alpha,gamma,beta) &
         &                +  qdrptens_cart(re,iatom,gamma,beta,alpha) )

           piezoim=piezoim+( qdrptens_cart(im,iatom,beta,alpha,gamma)-qdrptens_cart(im,iatom,alpha,gamma,beta) &
         &                +  qdrptens_cart(im,iatom,gamma,beta,alpha) )
         end if

       end do

       piezore=-piezore*half/ucvol
       piezoim=-piezoim*half/ucvol

       write(ab_out,'(3(i5,3x),2(1x,f20.10))') beta,gamma,alpha,piezore,piezoim

     end do
   end do
   write(ab_out,'(a)')' '
 end do
 write(ab_out,'(80a)')('=',ii=1,80)

 ABI_DEALLOCATE(qdrptens_cart)
 ABI_DEALLOCATE(dqpol_cart)
 ABI_DEALLOCATE(cartflg)

 DBG_EXIT("COLL")

end subroutine dfpt_qdrpout
!!***

!!****f* ABINIT/dfpt_flexo
!! NAME
!!  dfpt_flexo
!!
!! FUNCTION
!! This routine computes the elements of the flexoelectric tensor as: 
!!     --> Electronic contribution: second q-gradient of the second mixed derivative 
!!         of the energy w.r.t an electric field and a metric perturbation.
!!
!! COPYRIGHT
!!  Copyright (C) 2018 ABINIT group (MR,MS)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  codvsn=code version
!!  doccde(mband*nkpt*nsppol)=derivative of occupancies wrt the energy
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  dyewdq(2,3,natom,3,natom,3)= First q-gradient of Ewald part of the dynamical matrix
!!  dyewdqdq(2,3,natom,3,natom,3,3)= Second q-gradient of Ewald part of the dynamical matrix 
!!             sumed over the second sublattice
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  kxc(nfft,nkxc)=exchange and correlation kernel
!!  mpert=maximum number of perturbations for output processing
!!  mpi_enreg=information about MPI parallelization
!!  nattyp(ntypat)= # atoms of each type.
!!  nfft=(effective) number of FFT grid points (for this proc)
!!  ngfft(1:18)=integer array with FFT box dimensions and other
!!  nkpt=number of k points in the full BZ
!!  nkxc=second dimension of the kxc array. If /=0, the XC kernel must be computed.
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(mband*nkpt*nsppol)=occup number for each band (often 2) at each k point
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data for the GS (DUMMY)
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data (DUMMY)
!!  pertsy(3,natom+6)=set of perturbations that form a basis for all other perturbations
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  rhog(2,nfftf)=array for Fourier transform of GS electron density
!!  rhor(nfftf,nspden)=array for GS electron density in electrons/bohr**3.
!!  timrev=1 if time-reversal preserves the q wavevector; 0 otherwise.
!!  ucvol=unit cell volume in bohr**3.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  blkflg(3,mpert,3,mpert,3,mpert)= ( 1 if the element of the 3dte
!!   has been calculated ; 0 otherwise )
!!  d3etot(2,3,mpert,3,mpert,3,mpert)= matrix of the 3DTE
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!!  respfn
!!
!! CHILDREN
!!  appdig,dfpt_ciflexoout,dfpt_ciflexowf,dfpt_rhotov,distrb2,dotprod_vn,
!!  ebands_free,ebands_init,fourdp,
!!  getcut,getmpw,getph,hdr_init,hdr_update,
!!  initmpi_band,init_hamiltonian,initylmg,inwffil,kpgio
!!  load_spin_hamiltonian,hartredq,rf2_getidirs,read_rhor,
!!  setsym,symkpt,WffClose,wfk_open_read,wrtout,xmpi_sum
!!
!! SOURCE

subroutine dfpt_flexo(atindx,blkflg,codvsn,d3etot,doccde,dtfil,dtset,dyewdq,dyewdqdq,&
&          gmet,gprimd,kxc,mpert,&
&          mpi_enreg,nattyp,nfft,ngfft,nkpt,nkxc,&
&          nspden,nsppol,occ,pawrhoij,pawtab,pertsy,psps,rmet,rprimd,rhog,rhor,&
&          timrev,ucvol,xred)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_flexo'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpert,nfft,nkpt,nkxc,nspden,nsppol,timrev
 real(dp),intent(in) :: ucvol
 character(len=8), intent(in) :: codvsn
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(hdr_type) :: hdr_den
 type(pseudopotential_type),intent(in) :: psps
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)
!arrays
 integer,intent(in) :: atindx(dtset%natom)
 integer,intent(inout) :: blkflg(3,mpert,3,mpert,3,mpert)
 integer,intent(in) :: nattyp(dtset%ntypat),ngfft(18)
 integer,intent(in) :: pertsy(3,dtset%natom+6)
 real(dp),intent(inout) :: d3etot(2,3,mpert,3,mpert,3,mpert)
 real(dp),intent(in) :: doccde(dtset%mband*nkpt*dtset%nsppol)
 real(dp),intent(in) :: dyewdq(2,3,dtset%natom,3,dtset%natom,3)
 real(dp),intent(inout) :: dyewdqdq(2,3,dtset%natom,3,3,3)
 real(dp),intent(in) :: gmet(3,3), gprimd(3,3)
 real(dp),intent(in) :: kxc(nfft,nkxc)
 real(dp),intent(in) :: occ(dtset%mband*nkpt*dtset%nsppol)
 real(dp),intent(in) :: rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: rhog(2,nfft),rhor(nfft,nspden)
 real(dp),intent(in) :: xred(3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer :: ask_accurate,bantot_rbz,bdtot_index,bdtot1_index
 integer :: cplex,formeig,forunit,gscase,icg
 integer :: ia1,iatdir,iatom,iatpert,iatpert_cnt,iatpol
 integer :: ii,iefipert,iefipert_cnt,ierr,ikg,ikpt,ikpt1,ilm
 integer :: iq1dir,iq2dir,iq1grad,iq1grad_cnt,iq1q2grad,iq1q2grad_var
 integer :: ireadwf0,isppol,istrdir,istrpert,istrtype,istrpert_cnt,istwf_k,itypat,jatpert,jj,ka,kb
 integer :: lw_flexo,master,matom,matpert,mcg,me,mgfft,mkmem_rbz,mpw,my_nkpt_rbz
 integer :: natpert,nband_k,nefipert,nfftot,nhat1grdim,nkpt_rbz
 integer :: npw_k,npw1_k,nq1grad,nq1q2grad,nstrpert,nsym1,n3xccc
 integer :: nylmgr,optene,option,optorth,optres
 integer :: pawread,pertcase,qdir,spaceworld
 integer :: usexcnhat,useylmgr
 integer,parameter :: formeig1=1
 integer,parameter :: re=1,im=2
 real(dp) :: boxcut,delad,delag,delbd,delbg
 real(dp) :: doti,dotr,dum_scl,ecut_eff,ecut,etotal,fac,fermie,gsqcut,residm
 real(dp) :: vres2,wtk_k
 logical :: non_magnetic_xc
 character(len=500) :: msg                   
 character(len=fnlen) :: fi1o,fiwfatdis,fiwfstrain,fiwfefield,fiwfddk,fiwfdkdk
 type(ebands_t) :: bs_rbz
 type(hdr_type) :: hdr0
 type(wvl_data) :: wvl 
 type(wffile_type) :: wffgs,wfftgs
 type(gs_hamiltonian_type) :: gs_hamkq

!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 integer,allocatable :: bz2ibz_smap(:,:)
 integer,allocatable :: ddmdq_flg(:,:,:,:,:),elflexoflg(:,:,:,:)
 integer,allocatable :: indkpt1(:), indkpt1_tmp(:)
 integer,allocatable :: indsy1(:,:,:),irrzon1(:,:,:)
 integer,allocatable :: istwfk_rbz(:),isdq_flg(:,:,:,:,:)
 integer,allocatable :: kg(:,:),kg_k(:,:)
 integer,allocatable :: nband_rbz(:),npwarr(:),npwtot(:)
 integer,allocatable :: pert_atdis(:,:),pert_atdis_tmp(:,:)
 integer,allocatable :: pert_efield(:,:),pert_efield_tmp(:,:)
 integer,allocatable :: pert_strain(:,:),pert_strain_tmp(:,:)
 integer,allocatable :: q1grad(:,:),q1grad_tmp(:,:),q1q2grad(:,:)
 integer,allocatable :: symaf1(:),symrc1(:,:,:),symrl1(:,:,:)
 real(dp) :: kpoint(3)
 real(dp),allocatable :: cg(:,:),doccde_rbz(:)
 real(dp),allocatable :: ddmdq_qgradhart(:,:,:,:)
 real(dp),allocatable :: ddmdqwf(:,:,:,:),ddmdqwf_k(:,:,:,:)
 real(dp),allocatable :: ddmdqwf_t1(:,:,:,:),ddmdqwf_t1_k(:,:,:,:)
 real(dp),allocatable :: ddmdqwf_t2(:,:,:,:),ddmdqwf_t2_k(:,:,:,:)
 real(dp),allocatable :: ddmdqwf_t3(:,:,:,:),ddmdqwf_t3_k(:,:,:,:)
 real(dp),allocatable :: eigen0(:)
 real(dp),allocatable :: elflexowf(:,:,:,:,:),elflexowf_k(:,:,:,:,:)
 real(dp),allocatable :: elflexowf_t1(:,:,:,:,:),elflexowf_t1_k(:,:,:,:,:)
 real(dp),allocatable :: elflexowf_t2(:,:,:,:,:),elflexowf_t2_k(:,:,:,:,:)
 real(dp),allocatable :: elflexowf_t3(:,:,:,:,:),elflexowf_t3_k(:,:,:,:,:)
 real(dp),allocatable :: elflexowf_t4(:,:,:,:,:),elflexowf_t4_k(:,:,:,:,:)
 real(dp),allocatable :: elflexowf_t5(:,:,:,:,:),elflexowf_t5_k(:,:,:,:,:)
 real(dp),allocatable :: elqgradhart(:,:,:,:,:)
 real(dp),allocatable :: frwfdq(:,:,:,:,:,:),frwfdq_k(:,:,:,:,:,:)
 real(dp),allocatable :: isdq_qgradhart(:,:,:,:,:,:)
 real(dp),allocatable :: isdqwf(:,:,:,:,:,:),isdqwf_k(:,:,:,:,:,:)
 real(dp),allocatable :: isdqwf_t1(:,:,:,:,:,:),isdqwf_t1_k(:,:,:,:,:,:)
 real(dp),allocatable :: isdqwf_t2(:,:,:,:,:,:),isdqwf_t2_k(:,:,:,:,:,:)
 real(dp),allocatable :: isdqwf_t3(:,:,:,:,:,:),isdqwf_t3_k(:,:,:,:,:,:)
 real(dp),allocatable :: isdqwf_t4(:,:,:,:,:,:),isdqwf_t4_k(:,:,:,:,:,:)
 real(dp),allocatable :: isdqwf_t5(:,:,:,:,:,:),isdqwf_t5_k(:,:,:,:,:,:)
 real(dp),allocatable :: kpt_rbz(:,:)
 real(dp),allocatable :: nhat(:,:),nhat1(:,:),nhat1gr(:,:,:) 
 real(dp),allocatable :: occ_k(:),occ_rbz(:)
 real(dp),allocatable :: ph1d(:,:),phnons1(:,:,:) 
 real(dp),allocatable :: rhog1_tmp(:,:),rhog1_atdis(:,:,:)
 real(dp),allocatable :: rhog1_efield(:,:,:),rhor1_atdis(:,:,:),rhor1_tmp(:,:),rhor1_real(:,:)
 real(dp),allocatable :: rhor1_strain(:,:,:)
 real(dp),allocatable :: vhartr1(:),vhxc1_atdis(:,:),vhxc1_efield(:,:),vhxc1_strain(:,:)
 real(dp),allocatable :: vpsp1(:),vqgradhart(:),vresid1(:,:),vxc1(:,:)
 real(dp),allocatable :: dum_vxc(:,:)
 real(dp),allocatable :: wtk_folded(:), wtk_rbz(:)
 real(dp),allocatable,target :: vtrial1(:,:)
 real(dp),allocatable :: tnons1(:,:)
 real(dp),allocatable :: xccc3d1(:)
 real(dp),allocatable :: ylm(:,:),ylm_k(:,:),ylmgr(:,:,:),ylmgr_k(:,:,:)
 type(pawrhoij_type),allocatable :: pawrhoij_read(:)
 type(wfk_t),allocatable :: wfk_t_atdis(:),wfk_t_efield(:),wfk_t_ddk(:)
 type(wfk_t),allocatable :: wfk_t_dkdk(:),wfk_t_strain(:,:)
! *************************************************************************

 DBG_ENTER("COLL")

!Anounce start of flexoelectric tensor calculation
 write(msg, '(a,80a,a,a,a)' ) ch10,('=',ii=1,80),ch10,&
&   ' ==> Compute Flexoelectric Tensor Related Magnitudes <== ',ch10
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

!Init parallelism
 spaceworld=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt
 master=0

!Get FFT grid(s) sizes (be careful !) See NOTES in the comments at the beginning of respfn.F90
 mgfft=dtset%mgfft
 ecut=dtset%ecut
 ecut_eff=ecut*(dtset%dilatmx)**2

!Compute large sphere cut-off gsqcut
 call getcut(boxcut,ecut,gmet,gsqcut,dtset%iboxcut,std_out,dtset%qptn,ngfft)

!Various initializations
 lw_flexo=dtset%lw_flexo
 cplex=2-timrev
 matom=dtset%natom
 usexcnhat=0
 n3xccc=0
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 non_magnetic_xc=.true.

!Generate the 1-dimensional phases
 ABI_ALLOCATE(ph1d,(2,3*(2*dtset%mgfft+1)*dtset%natom))
 call getph(atindx,dtset%natom,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),ph1d,xred)

!################# PERTURBATIONS AND q-GRADIENTS LABELLING ############################

!Determine which atomic displacement, electric field, strain and q-gradient directions 
!have to be evaluated taking into account the perturbation symmetries
 matpert=dtset%natom*3

!Atomic displacement
 if (lw_flexo==1.or.lw_flexo==3.or.lw_flexo==4) then
   ABI_ALLOCATE(pert_atdis_tmp,(3,matpert))
   pert_atdis_tmp=0
   iatpert_cnt=0
   do iatpol=1,matom 
     do iatdir=1,3
       if (pertsy(iatdir,iatpol)==1) then
         iatpert_cnt=iatpert_cnt+1
         pert_atdis_tmp(1,iatpert_cnt)=iatpol                !atom displaced
         pert_atdis_tmp(2,iatpert_cnt)=iatdir                !direction of displacement
         pert_atdis_tmp(3,iatpert_cnt)=(iatpol-1)*3+iatdir   !like pertcase in dfpt_loopert.f90
       end if
     end do
   end do
   natpert=iatpert_cnt
   ABI_ALLOCATE(pert_atdis,(3,natpert))
   do iatpert=1,natpert
     pert_atdis(:,iatpert)=pert_atdis_tmp(:,iatpert)
   end do
   ABI_DEALLOCATE(pert_atdis_tmp)
 end if

!Electric field
 if (lw_flexo==1.or.lw_flexo==2) then
   ABI_ALLOCATE(pert_efield_tmp,(3,3))
   pert_efield_tmp=0
   iefipert_cnt=0
   do iefipert=1,3
     if (pertsy(iefipert,dtset%natom+2)==1) then
       iefipert_cnt=iefipert_cnt+1
       pert_efield_tmp(1,iefipert_cnt)=dtset%natom+2              !electric field pert
       pert_efield_tmp(2,iefipert_cnt)=iefipert                     !electric field direction
       pert_efield_tmp(3,iefipert_cnt)=matpert+3+iefipert           !like pertcase in dfpt_loopert.f90
     end if
   end do
   nefipert=iefipert_cnt
   ABI_ALLOCATE(pert_efield,(3,nefipert))
   do iefipert=1,nefipert
     pert_efield(:,iefipert)=pert_efield_tmp(:,iefipert)
   end do
   ABI_DEALLOCATE(pert_efield_tmp)
 end if

!ddk
 !The q1grad is related with the response to the ddk
 ABI_ALLOCATE(q1grad_tmp,(3,3))
 q1grad_tmp=0
 iq1grad_cnt=0
 do iq1grad=1,3
   if (pertsy(iq1grad,dtset%natom+1)==1) then
     iq1grad_cnt=iq1grad_cnt+1
     q1grad_tmp(1,iq1grad_cnt)=dtset%natom+1                         !ddk perturbation
     q1grad_tmp(2,iq1grad_cnt)=iq1grad                               !ddk direction
     q1grad_tmp(3,iq1grad_cnt)=matpert+iq1grad                       !like pertcase in dfpt_loopert.f90
   end if
 end do
 nq1grad=iq1grad_cnt
 ABI_ALLOCATE(q1grad,(3,nq1grad))
 do iq1grad=1,nq1grad
   q1grad(:,iq1grad)=q1grad_tmp(:,iq1grad)
 end do
 ABI_DEALLOCATE(q1grad_tmp)

!d2_dkdk
 !For the evaluation of the 2nd order q-gradient, the 9 directios are activated because 
 !currently the program calculates by defect all the components of the d2_dkdk perturbation.
 !TODO: This will have to be modified in the future when ABINIT enables to calculate specific
 !components of the d2_dkdk
 if (lw_flexo==1.or.lw_flexo==2) then
   nq1q2grad=9
   ABI_ALLOCATE(q1q2grad,(4,nq1q2grad))
   do iq1q2grad=1,nq1q2grad
     call rf2_getidirs(iq1q2grad,iq1dir,iq2dir)
     if (iq1dir==iq2dir) then
       q1q2grad(1,iq1q2grad)=dtset%natom+10                    !d2_dkdk perturbation diagonal elements
     else
       q1q2grad(1,iq1q2grad)=dtset%natom+11                    !d2_dkdk perturbation off-diagonal elements
     end if
     q1q2grad(2,iq1q2grad)=iq1dir                              !dq1 direction 
     q1q2grad(3,iq1q2grad)=iq2dir                              !dq2 direction 
     iq1q2grad_var=iq1q2grad
     if (iq1q2grad>6) iq1q2grad_var=iq1q2grad-3                !Lower=Upper diagonal triangle matrix
     q1q2grad(4,iq1q2grad)=iq1q2grad_var+(dtset%natom+6)*3     !like pertcase in dfpt_loopert.f90
   end do
 end if

!Strain perturbation
 if (lw_flexo==1.or.lw_flexo==2.or.lw_flexo==4) then
   ABI_ALLOCATE(pert_strain_tmp,(6,9))
   pert_strain_tmp=0
   !tmp uniaxial components
   istrpert_cnt=0
   do istrdir=1,3
     if (pertsy(istrdir,dtset%natom+3)==1) then 
       istrpert_cnt=istrpert_cnt+1
       ka=idx(2*istrdir-1);kb=idx(2*istrdir)
       pert_strain_tmp(1,istrpert_cnt)=dtset%natom+3        !Uniaxial strain perturbation
       pert_strain_tmp(2,istrpert_cnt)=istrdir              !Uniaxial strain case 
       pert_strain_tmp(3,istrpert_cnt)=ka                   !Strain direction 1
       pert_strain_tmp(4,istrpert_cnt)=kb                   !Strain direction 2
       pert_strain_tmp(5,istrpert_cnt)=matpert+6+istrdir    !like pertcase in dfpt_loopert.f90
       pert_strain_tmp(6,istrpert_cnt)=istrdir              !Indexing for the second q-gradient of the metric Hamiltonian
     end if
   end do
   !tmp shear components
   do istrdir=1,3
     if (pertsy(istrdir,dtset%natom+4)==1) then 
       istrpert_cnt=istrpert_cnt+1
       ka=idx(2*(istrdir+3)-1);kb=idx(2*(istrdir+3))
       pert_strain_tmp(1,istrpert_cnt)=dtset%natom+4        !Shear strain perturbation
       pert_strain_tmp(2,istrpert_cnt)=istrdir              !Shear strain case
       pert_strain_tmp(3,istrpert_cnt)=ka                   !Strain direction 1
       pert_strain_tmp(4,istrpert_cnt)=kb                   !Strain direction 2
       pert_strain_tmp(5,istrpert_cnt)=matpert+9+istrdir    !like pertcase in dfpt_loopert.f90
       pert_strain_tmp(6,istrpert_cnt)=3+istrdir            !Indexing for the second q-gradient of the metric Hamiltonian
     end if
   end do
   do istrdir=1,3
     if (pertsy(istrdir,dtset%natom+4)==1) then 
       istrpert_cnt=istrpert_cnt+1
       ka=idx(2*(istrdir+3)-1);kb=idx(2*(istrdir+3))
       pert_strain_tmp(1,istrpert_cnt)=dtset%natom+4        !Shear strain perturbation
       pert_strain_tmp(2,istrpert_cnt)=istrdir              !Shear strain case
       pert_strain_tmp(3,istrpert_cnt)=kb                   !Strain direction 1
       pert_strain_tmp(4,istrpert_cnt)=ka                   !Strain direction 2
       pert_strain_tmp(5,istrpert_cnt)=matpert+9+istrdir    !like pertcase in dfpt_loopert.f90
       pert_strain_tmp(6,istrpert_cnt)=6+istrdir            !Indexing for the second q-gradient of the metric Hamiltonian
     end if
   end do
   nstrpert=istrpert_cnt
   ABI_ALLOCATE(pert_strain,(6,nstrpert))
   do istrpert=1,nstrpert
     pert_strain(:,istrpert)=pert_strain_tmp(:,istrpert)
   end do
   ABI_DEALLOCATE(pert_strain_tmp)
end if

!################# ELECTROSTATIC CONTRIBUTIONS  #######################################

!This is necessary to deactivate paw options in the dfpt_rhotov routine
 ABI_DATATYPE_ALLOCATE(pawrhoij_read,(0))
 pawread=0
 nhat1grdim=0
 ABI_ALLOCATE(nhat1gr,(0,0,0))
 ABI_ALLOCATE(nhat,(nfft,nspden))
 nhat=zero
 ABI_ALLOCATE(nhat1,(cplex*nfft,nspden))
 nhat1=zero

!Read the first order densities response from a disk file, calculates the FFT 
!(rhog1_tmp) and the first order Hartree and xc potentials(vhxc1_pert). 
!TODO: In the call to read_rhor there is a security option that compares with the header
!hdr. Not activated at this moment.
 ABI_ALLOCATE(rhog1_tmp,(2,nfft))
 ABI_ALLOCATE(rhor1_tmp,(cplex*nfft,nspden))
 ABI_ALLOCATE(rhor1_real,(1*nfft,nspden))
 ABI_ALLOCATE(vhartr1,(cplex*nfft))
 ABI_ALLOCATE(vpsp1,(cplex*nfft))
 ABI_ALLOCATE(vtrial1,(cplex*nfft,nspden))
 ABI_ALLOCATE(vresid1,(cplex*nfft,nspden))
 ABI_ALLOCATE(vxc1,(cplex*nfft,nspden))
 ABI_ALLOCATE(dum_vxc,(nfft,nspden))
 ABI_ALLOCATE(xccc3d1,(cplex*n3xccc))
 vpsp1=zero; dum_vxc=zero
 optene=0; optres=1 

!Atomic displacement
 if (lw_flexo==1.or.lw_flexo==3.or.lw_flexo==4) then
   ABI_ALLOCATE(rhor1_atdis,(natpert,cplex*nfft,nspden))
   ABI_ALLOCATE(rhog1_atdis,(natpert,2,nfft))
   ABI_ALLOCATE(vhxc1_atdis,(natpert,cplex*nfft))
   vtrial1=zero
   do iatpert= 1, natpert
     iatpol=pert_atdis(1,iatpert)
     iatdir=pert_atdis(2,iatpert)
     pertcase=pert_atdis(3,iatpert)

     !Reads a real first order density
     call appdig(pertcase,dtfil%fildens1in,fi1o)
     call read_rhor(fi1o, 1, nspden, nfft, ngfft, pawread, mpi_enreg, rhor1_real, &
      & hdr_den, pawrhoij_read, spaceworld)

     !Perform FFT rhor1 to rhog1
     call fourdp(cplex,rhog1_tmp,rhor1_real,-1,mpi_enreg,nfft,1,ngfft,0)

     !Accumulate density in meaningful complex arrays
     if (timrev==0) then
       do ii=1,nfft
         jj=ii*2
         rhor1_tmp(jj-1,:)=rhor1_real(ii,:)
       end do
     else if (timrev==1) then
       rhor1_tmp(:,:)=rhor1_real(:,:)
     end if
     rhog1_atdis(iatpert,:,:)=rhog1_tmp(:,:)
     rhor1_atdis(iatpert,:,:)=rhor1_tmp(:,:)

     !Calculate first order Hartree and xc potentials
     call dfpt_rhotov(cplex,dum_scl,dum_scl,dum_scl,dum_scl,dum_scl, &
      & gsqcut,iatdir,iatpol,&
      & dtset%ixc,kxc,mpi_enreg,dtset%natom,nfft,ngfft,nhat,nhat1,nhat1gr,nhat1grdim,nkxc,&
      & nspden,n3xccc,non_magnetic_xc,optene,optres,dtset%qptn,rhog,rhog1_tmp,rhor,rhor1_tmp,&
      & rprimd,ucvol,psps%usepaw,usexcnhat,vhartr1,vpsp1,vresid1,vres2,vtrial1,dum_vxc,vxc1,&
      & xccc3d1,dtset%ixcrot)

     !Accumulate the potential in meaningful arrays
     vhxc1_atdis(iatpert,:)=vtrial1(:,nspden)

   end do
 end if

!Electric field
 if (lw_flexo==1.or.lw_flexo==2) then
   ABI_ALLOCATE(rhog1_efield,(nefipert,2,nfft))
   ABI_ALLOCATE(vhxc1_efield,(nefipert,cplex*nfft))
   vtrial1=zero
   do iefipert=1,nefipert
     pertcase=pert_efield(3,iefipert)
   
     !Reads a real first order density
     call appdig(pertcase,dtfil%fildens1in,fi1o)
     call read_rhor(fi1o, 1, nspden, nfft, ngfft, pawread, mpi_enreg, rhor1_real, &
      & hdr_den, pawrhoij_read, spaceworld)

     !Perform FFT rhor1 to rhog1
     call fourdp(cplex,rhog1_tmp,rhor1_real,-1,mpi_enreg,nfft,1,ngfft,0)

     !Accumulate density in meaningful complex arrays
     if (timrev==0) then
       do ii=1,nfft
         jj=ii*2
         rhor1_tmp(jj-1,:)=rhor1_real(ii,:)
       end do
     else if (timrev==1) then
       rhor1_tmp(:,:)=rhor1_real(:,:)
     end if
     rhog1_efield(iefipert,:,:)=rhog1_tmp(:,:)

     !Calculate first order Hartree and xc potentials
     call dfpt_rhotov(cplex,dum_scl,dum_scl,dum_scl,dum_scl,dum_scl, &
      & gsqcut,pert_efield(2,iefipert),dtset%natom+2,&
      & dtset%ixc,kxc,mpi_enreg,dtset%natom,nfft,ngfft,nhat,nhat1,nhat1gr,nhat1grdim,nkxc,&
      & nspden,n3xccc,non_magnetic_xc,optene,optres,dtset%qptn,rhog,rhog1_tmp,rhor,rhor1_tmp,&
      & rprimd,ucvol,psps%usepaw,usexcnhat,vhartr1,vpsp1,vresid1,vres2,vtrial1,dum_vxc,vxc1,&
      & xccc3d1,dtset%ixcrot)

     !Accumulate the potential in meaningful arrays
     vhxc1_efield(iefipert,:)=vtrial1(:,nspden)

   end do
endif

!Strain 
 if (lw_flexo==1.or.lw_flexo==2.or.lw_flexo==4) then
   ABI_ALLOCATE(rhor1_strain,(nstrpert,cplex*nfft,nspden))
   ABI_ALLOCATE(vhxc1_strain,(nstrpert,cplex*nfft))
   vtrial1=zero
   do istrpert= 1, nstrpert
     istrtype=pert_strain(1,istrpert)
     istrdir=pert_strain(2,istrpert)
     pertcase=pert_strain(5,istrpert)

     !Reads a real first order density
     call appdig(pertcase,dtfil%fildens1in,fi1o)
     call read_rhor(fi1o, 1, nspden, nfft, ngfft, pawread, mpi_enreg, rhor1_real, &
      & hdr_den, pawrhoij_read, spaceworld)

     !Perform FFT rhor1 to rhog1
     call fourdp(cplex,rhog1_tmp,rhor1_real,-1,mpi_enreg,nfft,1,ngfft,0)

     !Accumulate density in meaningful complex arrays
     if (timrev==0) then
       do ii=1,nfft
         jj=ii*2
         rhor1_tmp(jj-1,:)=rhor1_real(ii,:)
       end do
     else if (timrev==1) then
       rhor1_tmp(:,:)=rhor1_real(:,:)
     end if
     rhor1_strain(istrpert,:,:)=rhor1_tmp(:,:)

     !Calculate first order Hartree and xc potentials
     call dfpt_rhotov(cplex,dum_scl,dum_scl,dum_scl,dum_scl,dum_scl, &
      & gsqcut,istrdir,istrtype,&
      & dtset%ixc,kxc,mpi_enreg,dtset%natom,nfft,ngfft,nhat,nhat1,nhat1gr,nhat1grdim,nkxc,&
      & nspden,n3xccc,non_magnetic_xc,optene,optres,dtset%qptn,rhog,rhog1_tmp,rhor,rhor1_tmp,&
      & rprimd,ucvol,psps%usepaw,usexcnhat,vhartr1,vpsp1,vresid1,vres2,vtrial1,dum_vxc,vxc1,&
      & xccc3d1,dtset%ixcrot)

     !Accumulate the potential in meaningful arrays
     vhxc1_strain(istrpert,:)=vtrial1(:,nspden)

   end do
end if

 !These arrays will not be used anymore (for the moment)
 ABI_DEALLOCATE(rhor1_real)
 ABI_DEALLOCATE(rhor1_tmp)
 ABI_DEALLOCATE(vhartr1)
 ABI_DEALLOCATE(vpsp1)
 ABI_DEALLOCATE(vtrial1)
 ABI_DEALLOCATE(vresid1)
 ABI_DEALLOCATE(vxc1)
 ABI_DEALLOCATE(dum_vxc)
 ABI_DEALLOCATE(xccc3d1)

 ABI_DATATYPE_DEALLOCATE(pawrhoij_read)
 ABI_DEALLOCATE(nhat1gr)
 ABI_DEALLOCATE(nhat)
 ABI_DEALLOCATE(nhat1)

!!Calculate the electrostatic term from the q-gradient of the Hartree potential
 ABI_ALLOCATE(vqgradhart,(2*nfft))
 ABI_ALLOCATE(rhor1_tmp,(2*nfft,nspden))

!Electronic contribution 
 if (lw_flexo==1.or.lw_flexo==2) then
   ABI_ALLOCATE(elqgradhart,(2,3,3,3,3))
   ABI_ALLOCATE(elflexoflg,(3,3,3,3))
   elflexoflg=0
   rhor1_tmp=zero
   do iq1grad=1,nq1grad
     qdir=q1grad(2,iq1grad)
     do iefipert=1,nefipert

       !Calculate the gradient of the potential generated by the first order electric field density
       rhog1_tmp(:,:)=rhog1_efield(iefipert,:,:)
       call hartredq(2,gmet,gsqcut,mpi_enreg,nfft,ngfft,qdir,rhog1_tmp,vqgradhart) 

       !To ckeck
       !call appdig(pert_efield(3,iefipert)+q1grad(3,iq1grad),"Gradient_Hartree_potential",fi1o)
       !call fftdatar_write_from_hdr("first_order_potential",fi1o,dtset%iomode,hdr_den,&
       ! & ngfft,cplex,nfft,nspden,vqgradhart,mpi_enreg)

       do istrpert=1,nstrpert

         !Calculate the electrostatic energy term with the first order strain density 
         if (timrev==1) then
           do ii=1,nfft
             jj=ii*2
             rhor1_tmp(jj-1,:)=rhor1_strain(istrpert,ii,:)
           end do
         else if (timrev==0) then
           rhor1_tmp(:,:)=rhor1_strain(istrpert,:,:)
         end if
       
         call dotprod_vn(2,rhor1_tmp,dotr,doti,nfft,nfftot,nspden,2,vqgradhart,ucvol)
         elqgradhart(re,pert_efield(2,iefipert),q1grad(2,iq1grad),pert_strain(3,istrpert),pert_strain(4,istrpert))=dotr*half
         elqgradhart(im,pert_efield(2,iefipert),q1grad(2,iq1grad),pert_strain(3,istrpert),pert_strain(4,istrpert))=doti*half
         elflexoflg(pert_efield(2,iefipert),q1grad(2,iq1grad),pert_strain(3,istrpert),pert_strain(4,istrpert))=1

         blkflg(pert_efield(2,iefipert),pert_efield(1,iefipert), &
       &        pert_strain(2,istrpert),pert_strain(1,istrpert), &
       &        q1grad(2,iq1grad),matom+8)=1           

       end do
     end do
   end do 
 end if

!1st q-gradient of DM contribution
 if (lw_flexo==1.or.lw_flexo==3) then 
   ABI_ALLOCATE(ddmdq_qgradhart,(2,natpert,natpert,nq1grad))
   ABI_ALLOCATE(ddmdq_flg,(matom,3,matom,3,3))
   ddmdq_flg=0
   rhor1_tmp=zero
   do iq1grad=1,nq1grad
     qdir=q1grad(2,iq1grad)
     do iatpert=1,natpert

       rhog1_tmp(:,:)=rhog1_atdis(iatpert,:,:)
       call hartredq(2,gmet,gsqcut,mpi_enreg,nfft,ngfft,qdir,rhog1_tmp,vqgradhart) 

       !TODO:Maybe it is only necessary to compute half of these elements by symmetry
       do jatpert=1,natpert

         !Calculate the electrostatic energy term with the first order electric field density 
         if (timrev==1) then
           do ii=1,nfft
             jj=ii*2
             rhor1_tmp(jj-1,:)=rhor1_atdis(jatpert,ii,:)
           end do
         else if (timrev==0) then
           rhor1_tmp(:,:)=rhor1_atdis(jatpert,:,:)
         end if
       
         call dotprod_vn(2,rhor1_tmp,dotr,doti,nfft,nfftot,nspden,2,vqgradhart,ucvol)
         ddmdq_qgradhart(re,iatpert,jatpert,iq1grad)=dotr*half
         ddmdq_qgradhart(im,iatpert,jatpert,iq1grad)=doti*half
         ddmdq_flg(pert_atdis(1,iatpert),pert_atdis(2,iatpert),&
                 & pert_atdis(1,jatpert),pert_atdis(2,jatpert),q1grad(2,iq1grad))=1

         blkflg(pert_atdis(2,iatpert),pert_atdis(1,iatpert), &
       &        pert_atdis(2,jatpert),pert_atdis(1,jatpert), &
       &        q1grad(2,iq1grad),matom+8)=1           

       end do
     end do
   end do
 end if 

!1st g-gradient of internal strain tensor contribution
 if (lw_flexo==1.or.lw_flexo==4) then
   ABI_ALLOCATE(isdq_qgradhart,(2,matom,3,3,3,3))
   ABI_ALLOCATE(isdq_flg,(matom,3,3,3,3))
   isdq_flg=0
   rhor1_tmp=zero
   do iq1grad=1,nq1grad
     qdir=q1grad(2,iq1grad)
     do iatpert=1,natpert

       rhog1_tmp(:,:)=rhog1_atdis(iatpert,:,:)
       call hartredq(2,gmet,gsqcut,mpi_enreg,nfft,ngfft,qdir,rhog1_tmp,vqgradhart) 

       do istrpert=1,nstrpert

         !Calculate the electrostatic energy term with the first order strain density 
         if (timrev==1) then
           do ii=1,nfft
             jj=ii*2
             rhor1_tmp(jj-1,:)=rhor1_strain(istrpert,ii,:)
           end do
         else if (timrev==0) then
           rhor1_tmp(:,:)=rhor1_strain(istrpert,:,:)
         end if
       
         call dotprod_vn(2,rhor1_tmp,dotr,doti,nfft,nfftot,nspden,2,vqgradhart,ucvol)

         isdq_qgradhart(re,pert_atdis(1,iatpert),pert_atdis(2,iatpert),q1grad(2,iq1grad), &
       & pert_strain(3,istrpert),pert_strain(4,istrpert))=dotr*half
         isdq_qgradhart(im,pert_atdis(1,iatpert),pert_atdis(2,iatpert),q1grad(2,iq1grad), &
       & pert_strain(3,istrpert),pert_strain(4,istrpert))=doti*half
         isdq_flg(pert_atdis(1,iatpert),pert_atdis(2,iatpert),q1grad(2,iq1grad), &
       & pert_strain(3,istrpert),pert_strain(4,istrpert))=1

         blkflg(pert_atdis(2,iatpert),pert_atdis(1,iatpert), &
       &        pert_strain(2,istrpert),pert_strain(1,istrpert), &
       &        q1grad(2,iq1grad),matom+8)=1           
     
       end do
     end do
   end do
 end if

 ABI_DEALLOCATE(rhor1_tmp)
 ABI_DEALLOCATE(rhog1_tmp)
 ABI_DEALLOCATE(vqgradhart)
 if (lw_flexo==1.or.lw_flexo==3.or.lw_flexo==4) then 
   ABI_DEALLOCATE(rhog1_atdis)
   ABI_DEALLOCATE(rhor1_atdis)
 end if
 if (lw_flexo==1.or.lw_flexo==2) ABI_DEALLOCATE(rhog1_efield)
 if (lw_flexo==1.or.lw_flexo==2.or.lw_flexo==4) ABI_DEALLOCATE(rhor1_strain)

!################# WAVE FUNCTION CONTRIBUTIONS  #######################################

!Determine the subset of symmetry operations (nsym1 operations)
!that leaves the perturbation invariant, and initialize corresponding arrays
!symaf1, symrl1, tnons1 (and pawang1%zarot, if PAW)..
!MR TODO: For the moment only the identiy symmetry is activated (nsym1=1) 
!         In a future I will try to activate perturbation dependent symmetries
!         with littlegroup_pert.F90. 
 nsym1 = 1
 ABI_ALLOCATE(indsy1,(4,nsym1,dtset%natom))
 ABI_ALLOCATE(symrc1,(3,3,nsym1))
 ABI_ALLOCATE(symaf1,(nsym1))
 ABI_ALLOCATE(symrl1,(3,3,nsym1))
 ABI_ALLOCATE(tnons1,(3,nsym1))
 symaf1(1:nsym1)= 1
 symrl1(:,:,nsym1)= dtset%symrel(:,:,1)
 tnons1(:,nsym1)= 0_dp

!Set up corresponding symmetry data
 ABI_ALLOCATE(irrzon1,(dtset%nfft**(1-1/nsym1),2,(nspden/dtset%nsppol)-3*(nspden/4)))
 ABI_ALLOCATE(phnons1,(2,dtset%nfft**(1-1/nsym1),(nspden/dtset%nsppol)-3*(nspden/4)))
 call setsym(indsy1,irrzon1,1,dtset%natom,dtset%nfft,dtset%ngfft,nspden,dtset%nsppol,&
&nsym1,phnons1,symaf1,symrc1,symrl1,tnons1,dtset%typat,xred)

 ABI_DEALLOCATE(indsy1)
 ABI_DEALLOCATE(symaf1)
 ABI_DEALLOCATE(symrl1)
 ABI_DEALLOCATE(tnons1)

!Determine the subset of k-points needed in the "reduced Brillouin zone",
!and initialize other quantities
 ABI_ALLOCATE(indkpt1_tmp,(nkpt))
 ABI_ALLOCATE(wtk_folded,(nkpt))
 ABI_ALLOCATE(bz2ibz_smap,(6, nkpt))
 indkpt1_tmp(:)=0 

 if (dtset%kptopt==2) then
   call symkpt(0,gmet,indkpt1_tmp,ab_out,dtset%kptns,nkpt,nkpt_rbz,&
  & nsym1,symrc1,timrev,dtset%wtk,wtk_folded,bz2ibz_smap,xmpi_comm_self)
 else if (dtset%kptopt==3) then
   call symkpt(0,gmet,indkpt1_tmp,ab_out,dtset%kptns,nkpt,nkpt_rbz,&
  & nsym1,symrc1,0,dtset%wtk,wtk_folded,bz2ibz_smap,xmpi_comm_self)
 else
   write(msg,"(1a)") 'kptopt must be 2 or 3 for the quadrupole calculation'
   MSG_BUG(msg)
 end if
 ABI_DEALLOCATE(bz2ibz_smap)

 ABI_ALLOCATE(doccde_rbz,(dtset%mband*nkpt_rbz*dtset%nsppol))
 ABI_ALLOCATE(indkpt1,(nkpt_rbz))
 ABI_ALLOCATE(istwfk_rbz,(nkpt_rbz))
 ABI_ALLOCATE(kpt_rbz,(3,nkpt_rbz))
 ABI_ALLOCATE(nband_rbz,(nkpt_rbz*dtset%nsppol))
 ABI_ALLOCATE(occ_rbz,(dtset%mband*nkpt_rbz*dtset%nsppol))
 ABI_ALLOCATE(wtk_rbz,(nkpt_rbz))
 indkpt1(:)=indkpt1_tmp(1:nkpt_rbz)
 do ikpt=1,nkpt_rbz
     istwfk_rbz(ikpt)=dtset%istwfk(indkpt1(ikpt))
     kpt_rbz(:,ikpt)=dtset%kptns(:,indkpt1(ikpt))
     wtk_rbz(ikpt)=wtk_folded(indkpt1(ikpt))
 end do
 ABI_DEALLOCATE(indkpt1_tmp)
 ABI_DEALLOCATE(wtk_folded)

!Transfer occ to occ_rbz 
!NOTE : this takes into account that indkpt1 is ordered
!MG: What about using occ(band,kpt,spin) ???
 bdtot_index=0;bdtot1_index=0
 do isppol=1,dtset%nsppol
   ikpt1=1
   do ikpt=1,nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
!    Must test against ikpt1/=nkpt_rbz+1, before evaluate indkpt1(ikpt1)
     if(ikpt1/=nkpt_rbz+1)then
       if(ikpt==indkpt1(ikpt1))then
         nband_rbz(ikpt1+(isppol-1)*nkpt_rbz)=nband_k
         occ_rbz(1+bdtot1_index:nband_k+bdtot1_index)=occ(1+bdtot_index:nband_k+bdtot_index)
         doccde_rbz(1+bdtot1_index:nband_k+bdtot1_index)=doccde(1+bdtot_index:nband_k+bdtot_index)
         ikpt1=ikpt1+1
         bdtot1_index=bdtot1_index+nband_k
       end if
     end if
     bdtot_index=bdtot_index+nband_k
   end do
 end do

!Compute maximum number of planewaves at k
call getmpw(ecut_eff,dtset%exchn2n3d,gmet,istwfk_rbz,kpt_rbz,mpi_enreg,mpw,nkpt_rbz)

!Allocate some k-dependent arrays at k
 ABI_ALLOCATE(kg,(3,mpw*nkpt_rbz))
 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(npwarr,(nkpt_rbz))
 ABI_ALLOCATE(npwtot,(nkpt_rbz))

!Determine distribution of k-points/bands over MPI processes
 if (allocated(mpi_enreg%my_kpttab)) then
   ABI_DEALLOCATE(mpi_enreg%my_kpttab)
 end if
 ABI_ALLOCATE(mpi_enreg%my_kpttab,(nkpt_rbz))
 if(xmpi_paral==1) then
   ABI_ALLOCATE(mpi_enreg%proc_distrb,(nkpt_rbz,dtset%mband,dtset%nsppol))
   call distrb2(dtset%mband,nband_rbz,nkpt_rbz,mpi_enreg%nproc_cell,dtset%nsppol,mpi_enreg)
 else
   mpi_enreg%my_kpttab(:)=(/(ii,ii=1,nkpt_rbz)/)
 end if
 my_nkpt_rbz=maxval(mpi_enreg%my_kpttab)
 call initmpi_band(mpi_enreg,nband_rbz,nkpt_rbz,dtset%nsppol)
 mkmem_rbz =my_nkpt_rbz 
 
!Set up the basis sphere of planewaves at k
 call kpgio(ecut_eff,dtset%exchn2n3d,gmet,istwfk_rbz,kg,&
& kpt_rbz,mkmem_rbz,nband_rbz,nkpt_rbz,'PERS',mpi_enreg,mpw,npwarr,npwtot,dtset%nsppol)
 ABI_DEALLOCATE(npwtot)

!Set up the spherical harmonics (Ylm) and 1st gradients at k
 useylmgr=1; option=2 ; nylmgr=9
 ABI_ALLOCATE(ylm,(mpw*mkmem_rbz,psps%mpsang*psps%mpsang*psps%useylm))
 ABI_ALLOCATE(ylmgr,(mpw*mkmem_rbz,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr))
 if (psps%useylm==1) then
   call initylmg(gprimd,kg,kpt_rbz,mkmem_rbz,mpi_enreg,psps%mpsang,mpw,nband_rbz,nkpt_rbz,&
&   npwarr,dtset%nsppol,option,rprimd,ylm,ylmgr)
 end if

!Initialize band structure datatype at k
 bantot_rbz=sum(nband_rbz(1:nkpt_rbz*dtset%nsppol))
 ABI_ALLOCATE(eigen0,(bantot_rbz))
 eigen0(:)=zero
 call ebands_init(bantot_rbz,bs_rbz,dtset%nelect,doccde_rbz,eigen0,istwfk_rbz,kpt_rbz,&
& nband_rbz,nkpt_rbz,npwarr,dtset%nsppol,dtset%nspinor,dtset%tphysel,dtset%tsmear,dtset%occopt,occ_rbz,wtk_rbz,&
& dtset%charge, dtset%kptopt, dtset%kptrlatt_orig, dtset%nshiftk_orig, dtset%shiftk_orig, &
& dtset%kptrlatt, dtset%nshiftk, dtset%shiftk)
 ABI_DEALLOCATE(eigen0)
 ABI_DEALLOCATE(doccde_rbz)

!Initialize header, update it with evolving variables
 gscase=0 ! A GS WF file is read
 call hdr_init(bs_rbz,codvsn,dtset,hdr0,pawtab,gscase,psps,wvl%descr,&
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

 call hdr0%update(bantot_rbz,etotal,fermie,&
& residm,rprimd,occ_rbz,pawrhoij,xred,dtset%amu_orig(:,1),&
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

!Clean band structure datatype (should use it more in the future !)
 call ebands_free(bs_rbz)

!Initialize GS wavefunctions at k
 ireadwf0=1; formeig=0 ; ask_accurate=1 ; optorth=0
 mcg=mpw*dtset%nspinor*dtset%mband*mkmem_rbz*dtset%nsppol
 if (one*mpw*dtset%nspinor*dtset%mband*mkmem_rbz*dtset%nsppol > huge(1)) then
   write (msg,'(4a, 5(a,i0), 2a)')&
&   "Default integer is not wide enough to store the size of the GS wavefunction array (WF0, mcg).",ch10,&
&   "Action: increase the number of processors. Consider also OpenMP threads.",ch10,&
&   "nspinor: ",dtset%nspinor, "mpw: ",mpw, "mband: ",dtset%mband, "mkmem_rbz: ",&
&   mkmem_rbz, "nsppol: ",dtset%nsppol,ch10,&
&   'Note: Compiling with large int (int64) requires a full software stack (MPI/FFTW/BLAS/LAPACK...) compiled in int64 mode'
   MSG_ERROR(msg)
 end if
 ABI_STAT_ALLOCATE(cg,(2,mcg), ierr)
 ABI_CHECK(ierr==0, "out-of-memory in cg")

 ABI_ALLOCATE(eigen0,(dtset%mband*nkpt_rbz*dtset%nsppol))
 call inwffil(ask_accurate,cg,dtset,dtset%ecut,ecut_eff,eigen0,dtset%exchn2n3d,&
& formeig,hdr0,ireadwf0,istwfk_rbz,kg,&
& kpt_rbz,dtset%localrdwf,dtset%mband,mcg,&
& mkmem_rbz,mpi_enreg,mpw,nband_rbz,dtset%ngfft,nkpt_rbz,npwarr,&
& dtset%nsppol,dtset%nsym,occ_rbz,optorth,dtset%symafm,&
& dtset%symrel,dtset%tnons,dtfil%unkg,wffgs,wfftgs,&
& dtfil%unwffgs,dtfil%fnamewffk,wvl)
 ABI_DEALLOCATE(eigen0)
!Close wffgs%unwff, if it was ever opened (in inwffil)
 if (ireadwf0==1) then
   call WffClose(wffgs,ierr)
 end if

!==== Initialize most of the Hamiltonian ====
!1) Allocate all arrays and initialize quantities that do not depend on k and spin.
!2) Perform the setup needed for the non-local factors:
!3) Constant kleimann-Bylander energies are copied from psps to gs_hamkq.
 call init_hamiltonian(gs_hamkq,psps,pawtab,dtset%nspinor,nsppol,nspden,dtset%natom,&
& dtset%typat,xred,nfft,dtset%mgfft,ngfft,rprimd,dtset%nloalg,ph1d=ph1d,&
& use_gpu_cuda=dtset%use_gpu_cuda)


!==== Initialize response functions files and handlers ====
 !Atomic displacement
 if (lw_flexo==1.or.lw_flexo==3.or.lw_flexo==4) then
   ABI_ALLOCATE(wfk_t_atdis,(natpert))
   do iatpert=1,natpert
     pertcase=pert_atdis(3,iatpert)
     call appdig(pertcase,dtfil%fnamewff1,fiwfatdis)
  
     !The value 20 is taken arbitrarily I would say
     forunit=20+pertcase
  
     !Check that atdis file exists and open it
     if (.not. file_exists(fiwfatdis)) then
       ! Trick needed to run Abinit test suite in netcdf mode.
       if (file_exists(nctk_ncify(fiwfatdis))) then
         write(msg,"(3a)")"- File: ",trim(fiwfatdis),&
         " does not exist but found netcdf file with similar name."
         call wrtout(std_out,msg,'COLL')
         fiwfatdis = nctk_ncify(fiwfatdis)
       end if
       if (.not. file_exists(fiwfatdis)) then
         MSG_ERROR('Missing file: '//TRIM(fiwfatdis))
       end if
     end if
     write(msg,'(a,a)')'-open atomic displacement wf1 file :',trim(fiwfatdis)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
     call wfk_open_read(wfk_t_atdis(iatpert),fiwfatdis,formeig1,dtset%iomode,forunit,spaceworld)
     end do 
 end if

 !ddk files
 ABI_ALLOCATE(wfk_t_ddk,(nq1grad))
 do iq1grad=1,nq1grad
   pertcase=q1grad(3,iq1grad)
   call appdig(pertcase,dtfil%fnamewffddk,fiwfddk)

   !The value 20 is taken arbitrarily I would say
   forunit=20+pertcase

   !Check that ddk file exists and open it
   if (.not. file_exists(fiwfddk)) then
     ! Trick needed to run Abinit test suite in netcdf mode.
     if (file_exists(nctk_ncify(fiwfddk))) then
       write(msg,"(3a)")"- File: ",trim(fiwfddk),&
       " does not exist but found netcdf file with similar name."
       call wrtout(std_out,msg,'COLL')
       fiwfddk = nctk_ncify(fiwfddk)
     end if
     if (.not. file_exists(fiwfddk)) then
       MSG_ERROR('Missing file: '//TRIM(fiwfddk))
     end if
   end if
   write(msg, '(a,a)') '-open ddk wf1 file :',trim(fiwfddk)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   call wfk_open_read(wfk_t_ddk(iq1grad),fiwfddk,formeig1,dtset%iomode,forunit,spaceworld)
 end do 

 !Electric field files
 if (lw_flexo==1.or.lw_flexo==2) then
   ABI_ALLOCATE(wfk_t_efield,(nefipert))
   do iefipert=1,nefipert
     pertcase=pert_efield(3,iefipert)
     call appdig(pertcase,dtfil%fnamewff1,fiwfefield)

     !The value 20 is taken arbitrarily I would say
     forunit=20+pertcase

     !Check that efield file exists and open it
     if (.not. file_exists(fiwfefield)) then
       ! Trick needed to run Abinit test suite in netcdf mode.
       if (file_exists(nctk_ncify(fiwfefield))) then
         write(msg,"(3a)")"- File: ",trim(fiwfefield),&
         " does not exist but found netcdf file with similar name."
         call wrtout(std_out,msg,'COLL')
         fiwfefield = nctk_ncify(fiwfefield)
       end if
       if (.not. file_exists(fiwfefield)) then
         MSG_ERROR('Missing file: '//TRIM(fiwfefield))
       end if
     end if
     write(msg, '(a,a)') '-open electric field wf1 file :',trim(fiwfefield)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
     call wfk_open_read(wfk_t_efield(iefipert),fiwfefield,formeig1,dtset%iomode,forunit,spaceworld)
   end do
 end if

 !Strain files
 if (lw_flexo==1.or.lw_flexo==2.or.lw_flexo==4) then
   ABI_ALLOCATE(wfk_t_strain,(3,3))
   do istrpert=1,nstrpert
     pertcase=pert_strain(5,istrpert)
     ka=pert_strain(3,istrpert)
     kb=pert_strain(4,istrpert)
     call appdig(pertcase,dtfil%fnamewff1,fiwfstrain)
  
     !The value 20 is taken arbitrarily I would say
     forunit=20+pertcase
  
     !Check that strain file exists and open it
     if (.not. file_exists(fiwfstrain)) then
       ! Trick needed to run Abinit test suite in netcdf mode.
       if (file_exists(nctk_ncify(fiwfstrain))) then
         write(msg,"(3a)")"- File: ",trim(fiwfstrain),&
         " does not exist but found netcdf file with similar name."
         call wrtout(std_out,msg,'COLL')
         fiwfstrain = nctk_ncify(fiwfstrain)
       end if
       if (.not. file_exists(fiwfstrain)) then
         MSG_ERROR('Missing file: '//TRIM(fiwfstrain))
       end if
     end if
     if (ka>=kb) then
       write(msg, '(a,a)') '-open strain wf1 file :',trim(fiwfstrain)
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out,msg,'COLL')
       call wfk_open_read(wfk_t_strain(ka,kb),fiwfstrain,formeig1,dtset%iomode,forunit,spaceworld)
     else 
       wfk_t_strain(ka,kb)=wfk_t_strain(kb,ka)
     end if

   end do
 end if

 !d2_dkdk
 if (lw_flexo==1.or.lw_flexo==2) then
   ABI_ALLOCATE(wfk_t_dkdk,(nq1q2grad))
   do iq1q2grad=1,nq1q2grad
  
     pertcase=q1q2grad(4,iq1q2grad)
     call appdig(pertcase,dtfil%fnamewffdkdk,fiwfdkdk)
  
     !The value 20 is taken arbitrarily I would say
     forunit=20+pertcase
  
     !Check that d2_ddk file exists and open it
     if (.not. file_exists(fiwfdkdk)) then
       ! Trick needed to run Abinit test suite in netcdf mode.
       if (file_exists(nctk_ncify(fiwfdkdk))) then
         write(msg,"(3a)")"- File: ",trim(fiwfdkdk),&
         " does not exist but found netcdf file with similar name."
         call wrtout(std_out,msg,'COLL')
         fiwfdkdk = nctk_ncify(fiwfdkdk)
       end if
       if (.not. file_exists(fiwfdkdk)) then
         MSG_ERROR('Missing file: '//TRIM(fiwfdkdk))
       end if
     end if
     if (iq1q2grad <= 6) then
       write(msg, '(a,a)') '-open d2_dkdk wf2 file :',trim(fiwfdkdk)
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out,msg,'COLL')
       call wfk_open_read(wfk_t_dkdk(iq1q2grad),fiwfdkdk,formeig1,dtset%iomode,forunit,spaceworld)
     else
       wfk_t_dkdk(iq1q2grad)=wfk_t_dkdk(iq1q2grad-3)
     end if
   end do
 end if

!Allocate the electronic flexoelectric tensor part depending on the wave functions
 if (lw_flexo==1.or.lw_flexo==2) then
   ABI_ALLOCATE(elflexowf,(2,3,3,3,3))
   ABI_ALLOCATE(elflexowf_k,(2,3,3,3,3))
   ABI_ALLOCATE(elflexowf_t1,(2,3,3,3,3))
   ABI_ALLOCATE(elflexowf_t1_k,(2,3,3,3,3))
   ABI_ALLOCATE(elflexowf_t2,(2,3,3,3,3))
   ABI_ALLOCATE(elflexowf_t2_k,(2,3,3,3,3))
   ABI_ALLOCATE(elflexowf_t3,(2,3,3,3,3))
   ABI_ALLOCATE(elflexowf_t3_k,(2,3,3,3,3))
   ABI_ALLOCATE(elflexowf_t4,(2,3,3,3,3))
   ABI_ALLOCATE(elflexowf_t4_k,(2,3,3,3,3))
   ABI_ALLOCATE(elflexowf_t5,(2,3,3,3,3))
   ABI_ALLOCATE(elflexowf_t5_k,(2,3,3,3,3))
   elflexowf=zero
   elflexowf_t1=zero
   elflexowf_t2=zero
   elflexowf_t3=zero
   elflexowf_t4=zero
   elflexowf_t5=zero
 end if

!Allocate arrays for wf contributions to the first q-gradient of the dynamical matrix
 if (lw_flexo==1.or.lw_flexo==3) then 
   ABI_ALLOCATE(ddmdqwf,(2,natpert,natpert,nq1grad))
   ABI_ALLOCATE(ddmdqwf_k,(2,natpert,natpert,nq1grad))
   ABI_ALLOCATE(ddmdqwf_t1,(2,natpert,natpert,nq1grad))
   ABI_ALLOCATE(ddmdqwf_t1_k,(2,natpert,natpert,nq1grad))
   ABI_ALLOCATE(ddmdqwf_t2,(2,natpert,natpert,nq1grad))
   ABI_ALLOCATE(ddmdqwf_t2_k,(2,natpert,natpert,nq1grad))
   ABI_ALLOCATE(ddmdqwf_t3,(2,natpert,natpert,nq1grad))
   ABI_ALLOCATE(ddmdqwf_t3_k,(2,natpert,natpert,nq1grad))
   ddmdqwf=zero
   ddmdqwf_t1=zero
   ddmdqwf_t2=zero
   ddmdqwf_t3=zero
 end if

!Allocate arrays for wf contributions to the first q-gradient of the internal strain tensor 
 if (lw_flexo==1.or.lw_flexo==4) then 
   ABI_ALLOCATE(frwfdq,(2,matom,3,3,3,nq1grad))
   ABI_ALLOCATE(frwfdq_k,(2,matom,3,3,3,nq1grad))
   ABI_ALLOCATE(isdqwf,(2,matom,3,nq1grad,3,3))
   ABI_ALLOCATE(isdqwf_k,(2,matom,3,nq1grad,3,3))
   ABI_ALLOCATE(isdqwf_t1,(2,matom,3,nq1grad,3,3))
   ABI_ALLOCATE(isdqwf_t1_k,(2,matom,3,nq1grad,3,3))
   ABI_ALLOCATE(isdqwf_t2,(2,matom,3,nq1grad,3,3))
   ABI_ALLOCATE(isdqwf_t2_k,(2,matom,3,nq1grad,3,3))
   ABI_ALLOCATE(isdqwf_t3,(2,matom,3,nq1grad,3,3))
   ABI_ALLOCATE(isdqwf_t3_k,(2,matom,3,nq1grad,3,3))
   ABI_ALLOCATE(isdqwf_t4,(2,matom,3,3,3,nq1grad))
   ABI_ALLOCATE(isdqwf_t4_k,(2,matom,3,3,3,nq1grad))
   ABI_ALLOCATE(isdqwf_t5,(2,matom,3,nq1grad,3,3))
   ABI_ALLOCATE(isdqwf_t5_k,(2,matom,3,nq1grad,3,3))
   frwfdq=zero
   isdqwf=zero
   isdqwf_t1=zero
   isdqwf_t2=zero
   isdqwf_t3=zero
   isdqwf_t4=zero
   isdqwf_t5=zero
 end if

!LOOP OVER SPINS
 bdtot_index=0
 icg=0
 do isppol=1,nsppol
   ikg=0

!  Continue to initialize the Hamiltonian
   call gs_hamkq%load_spin(isppol,with_nonlocal=.true.)

!  BIG FAT k POINT LOOP
   do ikpt=1,nkpt_rbz

     nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
     istwf_k=istwfk_rbz(ikpt)
     npw_k=npwarr(ikpt)
     npw1_k=npw_k

     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
       bdtot_index=bdtot_index+nband_k

       cycle ! Skip the rest of the k-point loop
     end if

     ABI_ALLOCATE(occ_k,(nband_k))
     ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
     ABI_ALLOCATE(ylmgr_k,(npw_k,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr))
     occ_k(:)=occ_rbz(1+bdtot_index:nband_k+bdtot_index)
     kpoint(:)=kpt_rbz(:,ikpt)
     wtk_k=wtk_rbz(ikpt)

!    Get plane-wave vectors and related data at k
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     if (psps%useylm==1) then
       do ilm=1,psps%mpsang*psps%mpsang
         ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       end do
       if (useylmgr==1) then
         do ilm=1,psps%mpsang*psps%mpsang
           do ii=1,nylmgr
             ylmgr_k(1:npw_k,ii,ilm)=ylmgr(1+ikg:npw_k+ikg,ii,ilm)
           end do
         end do
       end if
     end if

!    Compute the wf contributions to the electronic flexoelectric tensor
     if (lw_flexo==1.or.lw_flexo==2) then
       call dfpt_ciflexowf(cg,cplex,dtset,elflexowf_k,elflexowf_t1_k,elflexowf_t2_k, &
       &  elflexowf_t3_k,elflexowf_t4_k,elflexowf_t5_k, &
       &  gs_hamkq,gsqcut,icg,ikpt,indkpt1,isppol,istwf_k, &
       &  kg_k,kpoint,mkmem_rbz, &
       &  mpi_enreg,mpw,nattyp,nband_k,nefipert,nfft,ngfft,nkpt_rbz, &
       &  npw_k,nq1grad,nq1q2grad,nspden,nsppol,nstrpert,nylmgr,occ_k, &
       &  pert_efield,pert_strain,ph1d,psps,q1grad,q1q2grad,rhog,rmet,ucvol,useylmgr, &
       &  vhxc1_efield,vhxc1_strain,wfk_t_efield,wfk_t_ddk, &
       &  wfk_t_dkdk,wfk_t_strain,wtk_k,ylm_k,ylmgr_k)

!      Add the contribution from each k-point
       elflexowf=elflexowf + elflexowf_k
       elflexowf_t1=elflexowf_t1 + elflexowf_t1_k
       elflexowf_t2=elflexowf_t2 + elflexowf_t2_k
       elflexowf_t3=elflexowf_t3 + elflexowf_t3_k
       elflexowf_t4=elflexowf_t4 + elflexowf_t4_k
       elflexowf_t5=elflexowf_t5 + elflexowf_t5_k
     end if

!    Compute the wf contributions to the first q-gradient of the dynamical matrix
     if (lw_flexo==1.or.lw_flexo==3) then
       call dfpt_ddmdqwf(atindx,cg,cplex,ddmdqwf_k,ddmdqwf_t1_k,ddmdqwf_t2_k,ddmdqwf_t3_k,&
       &  dtset,gs_hamkq,gsqcut,icg,ikpt,indkpt1,isppol,istwf_k,kg_k,kpoint,mkmem_rbz,    &
       &  mpi_enreg,mpw,natpert,nattyp,nband_k,nfft,ngfft,nkpt_rbz,npw_k,nq1grad,nspden,  &
       &  nsppol,nylmgr,occ_k,pert_atdis,ph1d,psps,q1grad,rmet,ucvol,useylmgr, &
       &  vhxc1_atdis,wfk_t_atdis,wfk_t_ddk,wtk_k,xred,ylm_k,ylmgr_k)

!      Add the contribution from each k-point
       ddmdqwf=ddmdqwf + ddmdqwf_k
       ddmdqwf_t1=ddmdqwf_t1 + ddmdqwf_t1_k
       ddmdqwf_t2=ddmdqwf_t2 + ddmdqwf_t2_k
       ddmdqwf_t3=ddmdqwf_t3 + ddmdqwf_t3_k
     end if

!    Compute the wf contributions to the first q-gradient of the internal strain tensor
     if (lw_flexo==1.or.lw_flexo==4) then

!      First calculate the frozen wf contribution (notice the type-I indexing of
!      this term)
       call dfpt_isdqfr(atindx,cg,cplex,dtset,frwfdq_k,gs_hamkq,gsqcut,icg,ikpt,&
       &  isppol,istwf_k,kg_k,kpoint,mkmem_rbz,mpi_enreg,matom,mpw,natpert,nattyp,nband_k,nfft,&
       &  ngfft,nkpt_rbz,npw_k,nq1grad,nspden,nsppol,nstrpert,nylmgr,occ_k,pert_atdis,   &
       &  pert_strain,ph1d,psps,rmet,ucvol,useylmgr,wtk_k,xred,ylm_k,ylmgr_k)

!      Add the contribution from each k-point
       frwfdq=frwfdq + frwfdq_k

!      Now comute the 1st order wf contributions
       call dfpt_isdqwf(atindx,cg,cplex,dtset,gs_hamkq,gsqcut,icg,ikpt,indkpt1,isdqwf_k, &
       &  isdqwf_t1_k,isdqwf_t2_k,isdqwf_t3_k,isdqwf_t4_k,isdqwf_t5_k,isppol,istwf_k, &
       &  kg_k,kpoint,matom,mkmem_rbz,mpi_enreg,mpw,natpert,nattyp,nband_k,nfft,ngfft,nkpt_rbz, &
       &  npw_k,nq1grad,nspden,nsppol,nstrpert,nylmgr,occ_k, &
       &  pert_atdis,pert_strain,ph1d,psps,q1grad,rhog,rmet,ucvol,useylmgr, &
       &  vhxc1_atdis,vhxc1_strain,wfk_t_atdis,wfk_t_ddk, &
       &  wfk_t_strain,wtk_k,xred,ylm_k,ylmgr_k)
       
!      Add the contribution from each k-point
       isdqwf=isdqwf + isdqwf_k
       isdqwf_t1=isdqwf_t1 + isdqwf_t1_k
       isdqwf_t2=isdqwf_t2 + isdqwf_t2_k
       isdqwf_t3=isdqwf_t3 + isdqwf_t3_k
       isdqwf_t4=isdqwf_t4 + isdqwf_t4_k
       isdqwf_t5=isdqwf_t5 + isdqwf_t5_k

     end if

!    Keep track of total number of bands
     bdtot_index=bdtot_index+nband_k

!    Shift arrays memory
     if (mkmem_rbz/=0) then
       icg=icg+npw_k*dtset%nspinor*nband_k
       ikg=ikg+npw_k
     end if

     ABI_DEALLOCATE(occ_k)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylmgr_k)

   end do
!  END BIG FAT k POINT LOOP
 end do
!END LOOP OVER SPINS

!Close response function files
 if (lw_flexo==1.or.lw_flexo==3.or.lw_flexo==4) then
   do iatpert=1,natpert
     call wfk_t_atdis(iatpert)%close()
   end do
 end if
 if (lw_flexo==1.or.lw_flexo==2.or.lw_flexo==4) then
   do istrpert=1,nstrpert
     if (istrpert <= 6) then
       ka=pert_strain(3,istrpert)
       kb=pert_strain(4,istrpert)
       call wfk_t_strain(ka,kb)%close()
     end if
   end do
 end if
 do iq1grad=1,nq1grad
   call wfk_t_ddk(iq1grad)%close()
 end do
 if (lw_flexo==1.or.lw_flexo==2) then
   do iefipert=1,nefipert
     call wfk_t_efield(iefipert)%close()
   end do
   do iq1q2grad=1,nq1q2grad
     if (iq1q2grad <= 6) call wfk_t_dkdk(iq1q2grad)%close()
   end do
 end if

!=== MPI communications ==================
 if (xmpi_paral==1) then

   if (lw_flexo==1.or.lw_flexo==2) then
     call xmpi_sum(elflexowf,spaceworld,ierr)
     call xmpi_sum(elflexowf_t1,spaceworld,ierr)
     call xmpi_sum(elflexowf_t2,spaceworld,ierr)
     call xmpi_sum(elflexowf_t3,spaceworld,ierr)
     call xmpi_sum(elflexowf_t4,spaceworld,ierr)
     call xmpi_sum(elflexowf_t5,spaceworld,ierr)
   end if

   if (lw_flexo==1.or.lw_flexo==3) then
     call xmpi_sum(ddmdqwf,spaceworld,ierr)
     call xmpi_sum(ddmdqwf_t1,spaceworld,ierr)
     call xmpi_sum(ddmdqwf_t2,spaceworld,ierr)
     call xmpi_sum(ddmdqwf_t3,spaceworld,ierr)
   end if

   if (lw_flexo==1.or.lw_flexo==4) then
     call xmpi_sum(frwfdq,spaceworld,ierr)
     call xmpi_sum(isdqwf,spaceworld,ierr)
     call xmpi_sum(isdqwf_t1,spaceworld,ierr)
     call xmpi_sum(isdqwf_t2,spaceworld,ierr)
     call xmpi_sum(isdqwf_t3,spaceworld,ierr)
     call xmpi_sum(isdqwf_t4,spaceworld,ierr)
     call xmpi_sum(isdqwf_t5,spaceworld,ierr)
   end if

 end if

!Sum the G=0 contribution to the geometric term of the first 
!q-gradient of the internal !strain tensor
 if (lw_flexo==1.or.lw_flexo==4) then
   fac=pi*pi

   !LOOP OVER ATOMIC DISPLACEMENT PERTURBATIONS
   do iatpert=1,natpert
     iatom=pert_atdis(1,iatpert)
     iatdir=pert_atdis(2,iatpert)

     !Determination of the atom type
     ia1=0
     itypat=0
     do ii=1,dtset%ntypat
       ia1=ia1+nattyp(ii)
       if(atindx(iatom)<=ia1.and.itypat==0)itypat=ii
     end do

     !LOOP OVER STRAIN PERTURBATIONS
     do istrpert= 1, nstrpert
       ka=pert_strain(3,istrpert)
       kb=pert_strain(4,istrpert)
       delad=zero ; if (iatdir==kb) delad=one
       delbd=zero ; if (ka==kb) delbd=one
  
       !LOOP OVER Q1-GRADIENT
       do iq1grad=1,3
         delag=zero ; if(iatdir==iq1grad) delag=one
         delbg=zero ; if(ka==iq1grad) delbg=one

         frwfdq(1,iatom,iatdir,ka,kb,iq1grad)=frwfdq(1,iatom,iatdir,ka,kb,iq1grad)+&
       & fac*rhog(1,1)*psps%vlspl(1,2,itypat)*(delag*delbd+delad*delbg)

       end do
     end do
   end do
 end if

!Anounce finalization of calculations
 if (lw_flexo==1.or.lw_flexo==2) then
   write(msg, '(a,a,a)' ) ch10, &
   ' Frozen-ion flexoelectric tensor calculation completed ',ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if
 if (lw_flexo==1.or.lw_flexo==3) then
   write(msg, '(a,a,a)' ) ch10, &
   ' Dynamical matrix 1st q-gradient calculation completed ',ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if
 if (lw_flexo==1.or.lw_flexo==4) then
   write(msg, '(a,a,a)' ) ch10, &
   ' Piezoelectric force response 1st q-gradient calculation completed ',ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if

!Gather the different terms in the flexoelectric tensor and print them out
 if (me==0) then
   if (lw_flexo==1.or.lw_flexo==2) then
     call dfpt_ciflexoout(d3etot,elflexoflg,elflexowf,elflexowf_t1,elflexowf_t2, &
   & elflexowf_t3,elflexowf_t4,elflexowf_t5, &
   & elqgradhart,gprimd,dtset%kptopt,matom,mpert,nefipert, &
   & nstrpert,nq1grad,pert_efield,pert_strain,dtset%prtvol,q1grad,rprimd,ucvol)
   end if
   if (lw_flexo==1.or.lw_flexo==3) then
     call dfpt_ddmdqout(ddmdq_flg,ddmdq_qgradhart,ddmdqwf,ddmdqwf_t1,ddmdqwf_t2,ddmdqwf_t3,d3etot,&
   & dyewdq,gprimd,dtset%kptopt,matom,mpert,natpert,nq1grad,pert_atdis,dtset%prtvol,q1grad,rprimd)
   end if
   if (lw_flexo==1.or.lw_flexo==4) then
     call dfpt_isdqout(d3etot,dyewdqdq,frwfdq,gprimd,isdq_flg,isdq_qgradhart,isdqwf,isdqwf_t1,isdqwf_t2,&
   & isdqwf_t3,isdqwf_t4,isdqwf_t5,dtset%kptopt,matom,mpert,natpert, &
   & nstrpert,nq1grad,pert_atdis,pert_strain,dtset%prtvol,q1grad,rprimd,ucvol)
   end if
 end if

!Deallocattions
 if (lw_flexo==1.or.lw_flexo==3.or.lw_flexo==4) then
   ABI_DEALLOCATE(pert_atdis)
   ABI_DEALLOCATE(vhxc1_atdis)
   ABI_DEALLOCATE(wfk_t_atdis)
 end if
 if (lw_flexo==1.or.lw_flexo==2) then
   ABI_DEALLOCATE(pert_efield)
   ABI_DEALLOCATE(q1q2grad)
   ABI_DEALLOCATE(vhxc1_efield)
   ABI_DEALLOCATE(elqgradhart)
   ABI_DEALLOCATE(elflexoflg)
   ABI_DEALLOCATE(wfk_t_efield)
   ABI_DEALLOCATE(wfk_t_dkdk)
   ABI_DEALLOCATE(elflexowf)
   ABI_DEALLOCATE(elflexowf_k)
   ABI_DEALLOCATE(elflexowf_t1)
   ABI_DEALLOCATE(elflexowf_t1_k)
   ABI_DEALLOCATE(elflexowf_t2)
   ABI_DEALLOCATE(elflexowf_t2_k)
   ABI_DEALLOCATE(elflexowf_t3)
   ABI_DEALLOCATE(elflexowf_t3_k)
   ABI_DEALLOCATE(elflexowf_t4)
   ABI_DEALLOCATE(elflexowf_t4_k)
   ABI_DEALLOCATE(elflexowf_t5)
   ABI_DEALLOCATE(elflexowf_t5_k)
 end if
 if (lw_flexo==1.or.lw_flexo==3) then 
   ABI_DEALLOCATE(ddmdq_qgradhart)
   ABI_DEALLOCATE(ddmdq_flg)
   ABI_DEALLOCATE(ddmdqwf)
   ABI_DEALLOCATE(ddmdqwf_k)
   ABI_DEALLOCATE(ddmdqwf_t1)
   ABI_DEALLOCATE(ddmdqwf_t1_k)
   ABI_DEALLOCATE(ddmdqwf_t2)
   ABI_DEALLOCATE(ddmdqwf_t2_k)
   ABI_DEALLOCATE(ddmdqwf_t3)
   ABI_DEALLOCATE(ddmdqwf_t3_k)
 end if 
 if (lw_flexo==1.or.lw_flexo==2.or.lw_flexo==4) then
   ABI_DEALLOCATE(pert_strain)
   ABI_DEALLOCATE(vhxc1_strain)
   ABI_DEALLOCATE(wfk_t_strain)
 end if
 if (lw_flexo==1.or.lw_flexo==4) then 
   ABI_DEALLOCATE(frwfdq)
   ABI_DEALLOCATE(frwfdq_k)
   ABI_DEALLOCATE(isdqwf)
   ABI_DEALLOCATE(isdqwf_k)
   ABI_DEALLOCATE(isdqwf_t1)
   ABI_DEALLOCATE(isdqwf_t1_k)
   ABI_DEALLOCATE(isdqwf_t2)
   ABI_DEALLOCATE(isdqwf_t2_k)
   ABI_DEALLOCATE(isdqwf_t3)
   ABI_DEALLOCATE(isdqwf_t3_k)
   ABI_DEALLOCATE(isdqwf_t4)
   ABI_DEALLOCATE(isdqwf_t4_k)
   ABI_DEALLOCATE(isdqwf_t5)
   ABI_DEALLOCATE(isdqwf_t5_k)
 end if
 ABI_DEALLOCATE(q1grad)
 ABI_DEALLOCATE(ph1d)
 ABI_DEALLOCATE(indkpt1)
 ABI_DEALLOCATE(istwfk_rbz)
 ABI_DEALLOCATE(kpt_rbz)
 ABI_DEALLOCATE(nband_rbz)
 ABI_DEALLOCATE(occ_rbz)
 ABI_DEALLOCATE(wtk_rbz)
 ABI_DEALLOCATE(kg)
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(npwarr)
 !ABI_DEALLOCATE(mpi_enreg%my_kpttab)
 !ABI_DEALLOCATE(mpi_enreg%proc_distrb)
 ABI_DEALLOCATE(ylm)
 ABI_DEALLOCATE(ylmgr)
 ABI_DEALLOCATE(wfk_t_ddk)
 if(xmpi_paral==1) then
   ABI_DEALLOCATE(mpi_enreg%proc_distrb)
 end if

 DBG_EXIT("COLL")

end subroutine dfpt_flexo
!!***

!!****f* ABINIT/dfpt_ciflexoout
!! NAME
!!  dfpt_ciflexoout
!!
!! FUNCTION
!!  This subroutine gathers the different terms entering the electrocic flexoelectric tensor,
!!  perfofms the transformation from reduced to cartesian coordinates and 
!!  writes out the tensor in output files.
!!  
!! COPYRIGHT
!!  Copyright (C) 2018 ABINIT group (MR,MS)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  elqgradhart(2,3,3,3,3)=electronic electrostatic contribution from the 
!!                                             q-gradient of the Hartree potential
!!  elflexoflg(3,3,3,3)=array that indicates which elements of the electronic contribution to the
!!                                             flexoelectric tensor have been calculated
!!  elflexowf(2,3,3,3,3)=total wave function contributions to the electronic flexoelectric tensor 
!!                       (except t4)
!!  elflexowf_t1(2,3,3,3,3)=term 1 of the wave function contribution 
!!  elflexowf_t2(2,3,3,3,3)=term 2 of the wave function contribution 
!!  elflexowf_t3(2,3,3,3,3)=term 3 of the wave function contribution 
!!  elflexowf_t4(2,3,3,3,3)=term 4 of the wave function contribution (type-I)
!!  elflexowf_t5(2,3,3,3,3)=term 5 of the wave function contribution 
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  kptopt=2 time reversal symmetry is enforced, 3 trs is not enforced (for debugging purposes)
!!  matom=number of atoms 
!!  mpert=maximum number of perturbations
!!  nefipert=number of electric field perturbations
!!  nstrpert=number of strain perturbations
!!  nq1grad=number of q1 (q_{\gamma}) gradients
!!  pert_efield(3,nefipert)=array with the info for the electric field perturbations
!!  pert_strain(6,nstrpert)=array with the info for the strain perturbations
!!  prtvol=volume of information to be printed. 1-> The different contributions to the quadrupole are printed.
!!  q1grad(3,nq1grad)=array with the info for the q1 (q_{\gamma}) gradients
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  ucvol=unit cell volume in bohr**3.
!!  
!! OUTPUT
!!  d3etot(2,3,mpert,3,mpert,3,mpert)= matrix of the 3DTE
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!!  dfpt_flexo
!!
!! CHILDREN
!!
!!  cart39 
!!
!! SOURCE

 subroutine dfpt_ciflexoout(d3etot,elflexoflg,elflexowf,elflexowf_t1,elflexowf_t2, &
    & elflexowf_t3,elflexowf_t4,elflexowf_t5, &
    & elqgradhart,gprimd,kptopt,matom,mpert,nefipert,&
    & nstrpert,nq1grad,pert_efield,pert_strain,prtvol,q1grad,rprimd,ucvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_ciflexoout'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: kptopt,matom,mpert,nefipert,nstrpert,nq1grad,prtvol
 real(dp),intent(inout) :: d3etot(2,3,mpert,3,mpert,3,mpert)
 real(dp),intent(in) :: ucvol

!arrays
 integer,intent(in) :: elflexoflg(3,3,3,3)
 integer,intent(in) :: pert_efield(3,nefipert)
 integer,intent(in) :: pert_strain(6,nstrpert)
 integer,intent(in) :: q1grad(3,nq1grad)
 real(dp),intent(in) :: elflexowf(2,3,3,3,3)
 real(dp),intent(inout) :: elflexowf_t1(2,3,3,3,3)
 real(dp),intent(inout) :: elflexowf_t2(2,3,3,3,3)
 real(dp),intent(inout) :: elflexowf_t3(2,3,3,3,3)
 real(dp),intent(inout) :: elflexowf_t4(2,3,3,3,3)
 real(dp),intent(inout) :: elflexowf_t5(2,3,3,3,3)
 real(dp),intent(inout) :: elqgradhart(2,3,3,3,3)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: rprimd(3,3)
 
!Local variables-------------------------------
!scalars
 integer :: alpha,beta,delta,efipert,gamma
 integer :: ibuf,iefidir,iefipert,ii,iq1dir,iq1grad,istr1dir,istr2dir,istrpert
 integer :: q1pert,strcomp,strpert
 integer, parameter :: re=1,im=2
 real(dp) :: fac,tmpim,tmpre,ucvolinv
 character(len=500) :: msg

!arrays
 integer,allocatable :: cartflg_t4(:,:,:,:),redflg(:,:,:,:)
 integer :: flg1(3),flg2(3)
 real(dp) :: vec1(3),vec2(3)
 real(dp),allocatable :: elec_flexotens_cart(:,:,:,:,:),elec_flexotens_red(:,:,:,:,:)
 real(dp),allocatable :: elflexowf_buffer_cart(:,:,:,:,:,:),elflexowf_t4_cart(:,:,:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

!Gather the different terms in the electronic contribution to the flexoelectric tensor
 ABI_ALLOCATE(elec_flexotens_red,(2,3,3,3,3))
! elec_flexotens_red=zero
 ucvolinv= 1.0_dp/ucvol

 if (kptopt==3) then

   !Compute real and 'true' imaginary parts of flexoelectric tensor and independent terms
   !T4 term needs further treatment and it will be lately added to the cartesian coordinates
   !version of the flexoelectric tensor
   do istrpert=1,nstrpert
     istr1dir=pert_strain(3,istrpert)
     istr2dir=pert_strain(4,istrpert)
     do iq1grad=1,nq1grad
       iq1dir=q1grad(2,iq1grad)
       do iefipert=1,nefipert
         iefidir=pert_efield(2,iefipert)

         if (elflexoflg(iefidir,iq1dir,istr1dir,istr2dir)==1) then  

           elec_flexotens_red(re,iefidir,iq1dir,istr1dir,istr2dir)=2.0_dp*ucvolinv* &
         & ( elqgradhart(re,iefidir,iq1dir,istr1dir,istr2dir) +               &
         &   elflexowf(re,iefidir,iq1dir,istr1dir,istr2dir) ) 
           elec_flexotens_red(im,iefidir,iq1dir,istr1dir,istr2dir)=2.0_dp*ucvolinv* &
         & ( elqgradhart(im,iefidir,iq1dir,istr1dir,istr2dir) +               &
         &   elflexowf(im,iefidir,iq1dir,istr1dir,istr2dir) ) 

           !Multiply by the imaginary unit that has been factorized out
           tmpre=elec_flexotens_red(re,iefidir,iq1dir,istr1dir,istr2dir)
           tmpim=elec_flexotens_red(im,iefidir,iq1dir,istr1dir,istr2dir)
           elec_flexotens_red(re,iefidir,iq1dir,istr1dir,istr2dir)=-tmpim
           elec_flexotens_red(im,iefidir,iq1dir,istr1dir,istr2dir)=tmpre

           !Divide by the imaginary unit the T4 term 
           !(this is because a -i factor has been factorized out in all the individual contributions to this therm)
           tmpre=elflexowf_t4(re,iefidir,istr1dir,istr2dir,iq1dir)
           tmpim=elflexowf_t4(im,iefidir,istr1dir,istr2dir,iq1dir)
           elflexowf_t4(re,iefidir,istr1dir,istr2dir,iq1dir)=tmpim*2.0_dp*ucvolinv
           elflexowf_t4(im,iefidir,istr1dir,istr2dir,iq1dir)=-tmpre*2.0_dp*ucvolinv
          
           !Compute and save individual terms in mixed coordinates
           if (prtvol==1) then

             tmpre=elqgradhart(re,iefidir,iq1dir,istr1dir,istr2dir)
             tmpim=elqgradhart(im,iefidir,iq1dir,istr1dir,istr2dir)
             elqgradhart(re,iefidir,iq1dir,istr1dir,istr2dir)=-tmpim*2.0_dp*ucvolinv
             elqgradhart(im,iefidir,iq1dir,istr1dir,istr2dir)=tmpre*2.0_dp*ucvolinv

             tmpre=elflexowf_t1(re,iefidir,iq1dir,istr1dir,istr2dir)
             tmpim=elflexowf_t1(im,iefidir,iq1dir,istr1dir,istr2dir)
             elflexowf_t1(re,iefidir,iq1dir,istr1dir,istr2dir)=-tmpim*2.0_dp*ucvolinv
             elflexowf_t1(im,iefidir,iq1dir,istr1dir,istr2dir)=tmpre*2.0_dp*ucvolinv

             tmpre=elflexowf_t2(re,iefidir,iq1dir,istr1dir,istr2dir)
             tmpim=elflexowf_t2(im,iefidir,iq1dir,istr1dir,istr2dir)
             elflexowf_t2(re,iefidir,iq1dir,istr1dir,istr2dir)=-tmpim*2.0_dp*ucvolinv
             elflexowf_t2(im,iefidir,iq1dir,istr1dir,istr2dir)=tmpre*2.0_dp*ucvolinv

             tmpre=elflexowf_t3(re,iefidir,iq1dir,istr1dir,istr2dir)
             tmpim=elflexowf_t3(im,iefidir,iq1dir,istr1dir,istr2dir)
             elflexowf_t3(re,iefidir,iq1dir,istr1dir,istr2dir)=-tmpim*2.0_dp*ucvolinv
             elflexowf_t3(im,iefidir,iq1dir,istr1dir,istr2dir)=tmpre*2.0_dp*ucvolinv

             tmpre=elflexowf_t5(re,iefidir,iq1dir,istr1dir,istr2dir)
             tmpim=elflexowf_t5(im,iefidir,iq1dir,istr1dir,istr2dir)
             elflexowf_t5(re,iefidir,iq1dir,istr1dir,istr2dir)=-tmpim*2.0_dp*ucvolinv
             elflexowf_t5(im,iefidir,iq1dir,istr1dir,istr2dir)=tmpre*2.0_dp*ucvolinv

           end if

         end if

       end do
     end do
   end do

 else if (kptopt==2) then
   
   !Compute real part of flexoelectric tensor and independent terms
   !T4 term needs further treatment and it will be lately added to the cartesian coordinates
   !version of the flexoelectric tensor
   do istrpert=1,nstrpert
     istr1dir=pert_strain(3,istrpert)
     istr2dir=pert_strain(4,istrpert)
     do iq1grad=1,nq1grad
       iq1dir=q1grad(2,iq1grad)
       do iefipert=1,nefipert
         iefidir=pert_efield(2,iefipert)

         if (elflexoflg(iefidir,iq1dir,istr1dir,istr2dir)==1) then

           elec_flexotens_red(re,iefidir,iq1dir,istr1dir,istr2dir)=2.0_dp*ucvolinv* &
         & ( elqgradhart(re,iefidir,iq1dir,istr1dir,istr2dir) +               &
         &   elflexowf(re,iefidir,iq1dir,istr1dir,istr2dir) ) 
           elec_flexotens_red(im,iefidir,iq1dir,istr1dir,istr2dir)=2.0_dp*ucvolinv* &
         & ( elqgradhart(im,iefidir,iq1dir,istr1dir,istr2dir) +               &
         &   elflexowf(im,iefidir,iq1dir,istr1dir,istr2dir) ) 

           !Multiply by the imaginary unit that has been factorized out
           tmpim=elec_flexotens_red(im,iefidir,iq1dir,istr1dir,istr2dir)
           elec_flexotens_red(re,iefidir,iq1dir,istr1dir,istr2dir)=-tmpim
           elec_flexotens_red(im,iefidir,iq1dir,istr1dir,istr2dir)=0.0_dp

           !Divide by the imaginary unit the T4 term 
           !(this is because a -i factor has been factorized out in all the individual contributions to this therm)
           tmpim=elflexowf_t4(im,iefidir,istr1dir,istr2dir,iq1dir)
           elflexowf_t4(re,iefidir,istr1dir,istr2dir,iq1dir)=2.0_dp*tmpim*ucvolinv
           elflexowf_t4(im,iefidir,istr1dir,istr2dir,iq1dir)=0.0_dp

           !Compute and save individual terms in mixed coordinates
           if (prtvol==1) then

             tmpim=elqgradhart(im,iefidir,iq1dir,istr1dir,istr2dir)
             elqgradhart(re,iefidir,iq1dir,istr1dir,istr2dir)=-tmpim*2.0_dp*ucvolinv
             elqgradhart(im,iefidir,iq1dir,istr1dir,istr2dir)=0.0_dp

             tmpim=elflexowf_t1(im,iefidir,iq1dir,istr1dir,istr2dir)
             elflexowf_t1(re,iefidir,iq1dir,istr1dir,istr2dir)=-tmpim*2.0_dp*ucvolinv
             elflexowf_t1(im,iefidir,iq1dir,istr1dir,istr2dir)=0.0_dp

             tmpim=elflexowf_t2(im,iefidir,iq1dir,istr1dir,istr2dir)
             elflexowf_t2(re,iefidir,iq1dir,istr1dir,istr2dir)=-tmpim*2.0_dp*ucvolinv
             elflexowf_t2(im,iefidir,iq1dir,istr1dir,istr2dir)=0.0_dp

             tmpim=elflexowf_t3(im,iefidir,iq1dir,istr1dir,istr2dir)
             elflexowf_t3(re,iefidir,iq1dir,istr1dir,istr2dir)=-tmpim*2.0_dp*ucvolinv
             elflexowf_t3(im,iefidir,iq1dir,istr1dir,istr2dir)=0.0_dp

             tmpim=elflexowf_t5(im,iefidir,iq1dir,istr1dir,istr2dir)
             elflexowf_t5(re,iefidir,iq1dir,istr1dir,istr2dir)=-tmpim*2.0_dp*ucvolinv
             elflexowf_t5(im,iefidir,iq1dir,istr1dir,istr2dir)=0.0_dp

           end if

         end if

       end do
     end do
   end do

 else
   write(msg,"(1a)") 'kptopt must be 2 or 3 for long-wave DFPT calculations'
   MSG_BUG(msg)
 end if

!Transormation to complete cartesian coordinates the flexoelectric tensor
!and separately the T4 term 
 ABI_ALLOCATE(elec_flexotens_cart,(2,3,3,3,3))
 ABI_ALLOCATE(elflexowf_t4_cart,(2,3,3,3,3))
 ABI_ALLOCATE(cartflg_t4,(3,3,3,3))
 elec_flexotens_cart=elec_flexotens_red
 elflexowf_t4_cart=elflexowf_t4
 cartflg_t4=0

! ABI_DEALLOCATE(elec_flexotens_red)

 if (prtvol==1) then
   ABI_ALLOCATE(elflexowf_buffer_cart,(5,2,3,3,3,3))
   elflexowf_buffer_cart(1,:,:,:,:,:)=elflexowf_t1(:,:,:,:,:)
   elflexowf_buffer_cart(2,:,:,:,:,:)=elflexowf_t2(:,:,:,:,:)
   elflexowf_buffer_cart(3,:,:,:,:,:)=elflexowf_t3(:,:,:,:,:)
   elflexowf_buffer_cart(4,:,:,:,:,:)=elqgradhart(:,:,:,:,:)
   elflexowf_buffer_cart(5,:,:,:,:,:)=elflexowf_t5(:,:,:,:,:)
 end if

!1st transform coordinates of the electric field derivative of the flexoelectric tensor
 do istr2dir=1,3
   do istr1dir=1,3
     do iq1dir=1,3
       do ii=1,2
         do iefidir=1,3
           vec1(iefidir)=elec_flexotens_cart(ii,iefidir,iq1dir,istr1dir,istr2dir)
           flg1(iefidir)=elflexoflg(iefidir,iq1dir,istr1dir,istr2dir)
         end do
         call cart39(flg1,flg2,gprimd,matom+2,matom,rprimd,vec1,vec2)
         do iefidir=1,3
           elec_flexotens_cart(ii,iefidir,iq1dir,istr1dir,istr2dir)=vec2(iefidir)
         end do
       end do
     end do
   end do
 end do

!Do now the transformation of the electric field derivative for the T4 term
 do istr2dir=1,3
   do istr1dir=1,3
     do iq1dir=1,3
       do ii=1,2
         do iefidir=1,3
           vec1(iefidir)=elflexowf_t4_cart(ii,iefidir,istr1dir,istr2dir,iq1dir)
           flg1(iefidir)=elflexoflg(iefidir,iq1dir,istr1dir,istr2dir)
         end do
         call cart39(flg1,flg2,gprimd,matom+2,matom,rprimd,vec1,vec2)
         do iefidir=1,3
           elflexowf_t4_cart(ii,iefidir,istr1dir,istr2dir,iq1dir)=vec2(iefidir)
           cartflg_t4(iefidir,istr1dir,istr2dir,iq1dir)=flg2(iefidir)
         end do
       end do
     end do
   end do
 end do

!Do now the transformation of the electric field derivative for the other therms
 if (prtvol==1) then
   do ibuf=1,5
     do istr2dir=1,3
       do istr1dir=1,3
         do iq1dir=1,3
           do ii=1,2
             do iefidir=1,3
               vec1(iefidir)=elflexowf_buffer_cart(ibuf,ii,iefidir,iq1dir,istr1dir,istr2dir)
               flg1(iefidir)=elflexoflg(iefidir,iq1dir,istr1dir,istr2dir)
             end do
             call cart39(flg1,flg2,gprimd,matom+2,matom,rprimd,vec1,vec2)
             do iefidir=1,3
               elflexowf_buffer_cart(ibuf,ii,iefidir,iq1dir,istr1dir,istr2dir)=vec2(iefidir)
             end do
           end do
         end do
       end do
     end do
   end do
 end if

!2nd transform coordinates of the q-gradient (treat it as electric field)
!of the flexoelectric tensor
 do istr2dir=1,3
   do istr1dir=1,3
     do iefidir=1,3
       do ii=1,2
         do iq1dir=1,3
           vec1(iq1dir)=elec_flexotens_cart(ii,iefidir,iq1dir,istr1dir,istr2dir)
           flg1(iq1dir)=elflexoflg(iefidir,iq1dir,istr1dir,istr2dir)
         end do
         call cart39(flg1,flg2,gprimd,matom+2,matom,rprimd,vec1,vec2)
         do iq1dir=1,3
           elec_flexotens_cart(ii,iefidir,iq1dir,istr1dir,istr2dir)=vec2(iq1dir)
         end do
       end do
     end do
   end do
 end do

!2nd transform coordinates of the q-gradient (treat it as electric field)
!of the individual terms
 if (prtvol==1) then
   do ibuf=1,5
     do istr2dir=1,3
       do istr1dir=1,3
         do iefidir=1,3
           do ii=1,2
             do iq1dir=1,3
               vec1(iq1dir)=elflexowf_buffer_cart(ibuf,ii,iefidir,iq1dir,istr1dir,istr2dir)
               flg1(iq1dir)=elflexoflg(iefidir,iq1dir,istr1dir,istr2dir)
             end do
             call cart39(flg1,flg2,gprimd,matom+2,matom,rprimd,vec1,vec2)
             do iq1dir=1,3
               elflexowf_buffer_cart(ibuf,ii,iefidir,iq1dir,istr1dir,istr2dir)=vec2(iq1dir)
             end do
           end do
         end do
       end do
     end do
   end do
 end if

!Write the flexoelectric tensor in cartesian coordinates
!Open output files
 if (prtvol==1) then
   open(unit=71,file='elec_flexo_wf_t1.out',status='unknown',form='formatted',action='write')
   open(unit=72,file='elec_flexo_wf_t2.out',status='unknown',form='formatted',action='write')
   open(unit=73,file='elec_flexo_wf_t3.out',status='unknown',form='formatted',action='write')
   open(unit=74,file='elec_flexo_wf_t4.out',status='unknown',form='formatted',action='write')
   open(unit=75,file='elec_flexo_wf_t5.out',status='unknown',form='formatted',action='write')
   open(unit=76,file='elec_flexo_elecstic.out',status='unknown',form='formatted',action='write')
 end if

 write(ab_out,'(a)')' '
 write(ab_out,'(a)')' Frozen-ion flexoelectric tensor, in cartesian coordinates,'
 write(ab_out,'(a)')' efidir  qgrdir  strdir1  strdir2         real part          imaginary part'
 do istr2dir=1,3
   delta=istr2dir
   do istr1dir=1,3
     beta=istr1dir
     do iq1dir=1,3
       gamma=iq1dir
       do iefidir=1,3
         alpha=iefidir

         if (cartflg_t4(alpha,beta,delta,gamma)==1 .and. cartflg_t4(alpha,delta,gamma,beta)==1 &
         & .and. cartflg_t4(alpha,gamma,beta,delta)==1) then

           !Converts the T4 term to type-II form
           tmpre= elflexowf_t4_cart(re,alpha,beta,delta,gamma) + &
                & elflexowf_t4_cart(re,alpha,delta,gamma,beta) - &
                & elflexowf_t4_cart(re,alpha,gamma,beta,delta)
           
           tmpim= elflexowf_t4_cart(im,alpha,beta,delta,gamma) + &
                & elflexowf_t4_cart(im,alpha,delta,gamma,beta) - &
                & elflexowf_t4_cart(im,alpha,gamma,beta,delta)

           !Add the T4 term after conversion to type-II form
           elec_flexotens_cart(re,alpha,gamma,beta,delta)= &
         & elec_flexotens_cart(re,alpha,gamma,beta,delta) + tmpre
           elec_flexotens_cart(im,alpha,gamma,beta,delta)= &
         & elec_flexotens_cart(im,alpha,gamma,beta,delta) + tmpim
         
           !Writes the complete flexoelectric tensor
           write(ab_out,'(4(i5,3x),2(1x,f20.10))') alpha,gamma,beta,delta, &
         & elec_flexotens_cart(re,alpha,gamma,beta,delta), &
         & elec_flexotens_cart(im,alpha,gamma,beta,delta)

           if (prtvol==1) then
             write(71,'(4(i5,3x),2(1x,f20.10))') alpha,gamma,beta,delta, &
           & elflexowf_buffer_cart(1,re,alpha,gamma,beta,delta), &
           & elflexowf_buffer_cart(1,im,alpha,gamma,beta,delta)
 
             write(72,'(4(i5,3x),2(1x,f20.10))') alpha,gamma,beta,delta, &
           & elflexowf_buffer_cart(2,re,alpha,gamma,beta,delta), &
           & elflexowf_buffer_cart(2,im,alpha,gamma,beta,delta)
 
             write(73,'(4(i5,3x),2(1x,f20.10))') alpha,gamma,beta,delta, &
           & elflexowf_buffer_cart(3,re,alpha,gamma,beta,delta), &
           & elflexowf_buffer_cart(3,im,alpha,gamma,beta,delta)
 
             write(74,'(4(i5,3x),2(1x,f20.10))') alpha,gamma,beta,delta, &
           & tmpre, tmpim

             write(75,'(4(i5,3x),2(1x,f20.10))') alpha,gamma,beta,delta, &
           & elflexowf_buffer_cart(5,re,alpha,gamma,beta,delta), &
           & elflexowf_buffer_cart(5,im,alpha,gamma,beta,delta)
 
             write(76,'(4(i5,3x),2(1x,f20.10))') alpha,gamma,beta,delta, &
           & elflexowf_buffer_cart(4,re,alpha,gamma,beta,delta), &
           & elflexowf_buffer_cart(4,im,alpha,gamma,beta,delta)
           end if

         end if
          
       end do
     end do
     write(ab_out,'(a)')' '
     if (prtvol==1) then
       write(71,*)' ' 
       write(72,*)' ' 
       write(73,*)' ' 
       write(74,*)' ' 
       write(75,*)' ' 
       write(76,*)' ' 
     end if
   end do
 end do

 if (prtvol==1) then
   close(71)
   close(72)
   close(73)
   close(74)
   close(75)
   close(76)
   ABI_DEALLOCATE(elflexowf_buffer_cart)
 end if

!Calculate the contribution to the d3etot in mixed (reduced/cartesian) coordinates
 elec_flexotens_red=elec_flexotens_cart
 ABI_DEALLOCATE(elec_flexotens_cart)
 ABI_DEALLOCATE(elflexowf_t4_cart)
 ABI_DEALLOCATE(cartflg_t4)
 ABI_ALLOCATE(redflg,(3,3,3,3))
 redflg=0

!1st transform back coordinates of the electric field derivative of the flexoelectric tensor
 fac=two_pi ** 2
 do istr2dir=1,3
   do istr1dir=1,3
     do iq1dir=1,3
       do ii=1,2
         do iefidir=1,3
           vec1(iefidir)=elec_flexotens_red(ii,iefidir,iq1dir,istr1dir,istr2dir)
           flg1(iefidir)=elflexoflg(iefidir,iq1dir,istr1dir,istr2dir)
         end do
         call cart39(flg1,flg2,transpose(rprimd),matom+2,matom,transpose(gprimd),vec1,vec2)
         do iefidir=1,3
           elec_flexotens_red(ii,iefidir,iq1dir,istr1dir,istr2dir)=vec2(iefidir)*fac
           redflg(iefidir,iq1dir,istr1dir,istr2dir)=flg2(iefidir)
         end do
       end do
     end do
   end do
 end do

!2nd transform back coordinates of the q-gradient (treat it as electric field)
!of the flexoelectric tensor
 do istr2dir=1,3
   do istr1dir=1,3
     do iefidir=1,3
       do ii=1,2
         do iq1dir=1,3
           vec1(iq1dir)=elec_flexotens_red(ii,iefidir,iq1dir,istr1dir,istr2dir)
           flg1(iq1dir)=elflexoflg(iefidir,iq1dir,istr1dir,istr2dir)
         end do
         call cart39(flg1,flg2,transpose(rprimd),matom+2,matom,transpose(gprimd),vec1,vec2)
         do iq1dir=1,3
           elec_flexotens_red(ii,iefidir,iq1dir,istr1dir,istr2dir)=vec2(iq1dir)*fac
         end do
       end do
     end do
   end do
 end do

!Add contributions to d3etot
 efipert=matom+2
 q1pert=matom+8
 fac=ucvol/two
 do istrpert=1,nstrpert
   strpert=pert_strain(1,istrpert)
   strcomp=pert_strain(2,istrpert)
   istr1dir=pert_strain(3,istrpert)
   istr2dir=pert_strain(4,istrpert)
   do iq1grad=1,nq1grad
     iq1dir=q1grad(2,iq1grad)
     do iefipert=1,nefipert
       iefidir=pert_efield(2,iefipert)

       if (redflg(iefidir,iq1dir,istr1dir,istr2dir)==1) then
         d3etot(re,iefidir,efipert,strcomp,strpert,iq1dir,q1pert)= &
       & elec_flexotens_red(im,iefidir,iq1dir,istr1dir,istr2dir)*fac
         d3etot(im,iefidir,efipert,strcomp,strpert,iq1dir,q1pert)= &
       & -elec_flexotens_red(re,iefidir,iq1dir,istr1dir,istr2dir)*fac
       end if
       
     end do
   end do
 end do

 ABI_DEALLOCATE(elec_flexotens_red)
 ABI_DEALLOCATE(redflg)


 DBG_EXIT("COLL")

end subroutine dfpt_ciflexoout
!!***


!!****f* ABINIT/dfpt_ddmdqout
!! NAME
!!  dfpt_ddmdqout
!!
!! FUNCTION
!!  This subroutine gathers the different terms entering the first q derivative of the dynamical 
!!  matrix, perfofms the transformation from reduced to cartesian coordinates and 
!!  writes out the tensor in output files.
!!  
!! COPYRIGHT
!!  Copyright (C) 2018 ABINIT group (MR,MS)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ddmdq_flg(natpert,natpert,nq1grad)=array that indicates which elements of the first q derivative
!!                                             of dynamical matrix have been calculated
!!  ddmdq_qgradhart(2,natpert,natpert,nq1grad)=electronic electrostatic contribution from the 
!!                                             q-gradient of the Hartree potential
!!  ddmdqwf(2,natpert,natpert,nq1grad)=total wave function contribution to the first q derivative of dynamical matrix
!!  ddmdqwf_t1(2,natpert,natpert,nq1grad)=term 1 of the wave function contribution 
!!  ddmdqwf_t2(2,natpert,natpert,nq1grad)=term 2 of the wave function contribution 
!!  ddmdqwf_t3(2,natpert,natpert,nq1grad)=term 3 of the wave function contribution 
!!  dyewdq(2,3,natom,3,natom,3)= First q-gradient of Ewald part of the dynamical matrix
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  kptopt=2 time reversal symmetry is enforced, 3 trs is not enforced (for debugging purposes)
!!  matom=number of atoms 
!!  mpert=maximum number of perturbations
!!  natpert=number of atomic displacement perturbations
!!  nq1grad=number of q1 (q_{\gamma}) gradients
!!  pert_atdis(3,natpert)=array with the info for the electric field perturbations
!!  prtvol=volume of information to be printed. 1-> The different contributions to the quadrupole are printed.
!!  q1grad(3,nq1grad)=array with the info for the q1 (q_{\gamma}) gradients
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  
!! OUTPUT
!!  d3etot(2,3,mpert,3,mpert,3,mpert)= matrix of the 3DTE
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!!  dfpt_flexo
!!
!! CHILDREN
!!
!!  cart39 
!!
!! SOURCE

 subroutine dfpt_ddmdqout(ddmdq_flg,ddmdq_qgradhart,ddmdqwf,ddmdqwf_t1,ddmdqwf_t2,ddmdqwf_t3,d3etot, &
 & dyewdq,gprimd,kptopt,matom,mpert,natpert,nq1grad,pert_atdis,prtvol,q1grad,rprimd)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_ciflexoout'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: kptopt,matom,mpert,natpert,nq1grad,prtvol
 real(dp),intent(inout) :: d3etot(2,3,mpert,3,mpert,3,mpert)

!arrays
 integer,intent(in) :: ddmdq_flg(natpert,natpert,nq1grad)
 integer,intent(in) :: pert_atdis(3,natpert)
 integer,intent(in) :: q1grad(3,nq1grad)
 real(dp),intent(inout) :: ddmdq_qgradhart(2,natpert,natpert,nq1grad)
 real(dp),intent(in) :: ddmdqwf(2,natpert,natpert,nq1grad)
 real(dp),intent(inout) :: ddmdqwf_t1(2,natpert,natpert,nq1grad)
 real(dp),intent(inout) :: ddmdqwf_t2(2,natpert,natpert,nq1grad)
 real(dp),intent(inout) :: ddmdqwf_t3(2,natpert,natpert,nq1grad)
 real(dp),intent(in) :: dyewdq(2,3,matom,3,matom,3)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: rprimd(3,3)
 
!Local variables-------------------------------
!scalars
 integer :: iatdir,iatom,iatpert,ii,iq1dir,iq1grad,iq1pert,jatdir,jatom,jatpert
 integer, parameter :: im=2,re=1
 real(dp) :: piezofrim,piezofrre,tmpim,tmpre
 character(len=500) :: msg

!arrays
 integer :: flg1(3),flg2(3)
 integer, allocatable :: cartflg(:,:,:,:,:),ddmdq_cartflg(:,:,:,:,:)
 real(dp) :: vec1(3),vec2(3)
 real(dp), allocatable :: ddmdq_cart(:,:,:,:,:,:),ddmdq_red(:,:,:,:)

!****************************************************************************

 DBG_ENTER("COLL")

!Open output files
 if (prtvol==1) then
  open(unit=71,file='ddmdq_wf_t1.out',status='unknown',form='formatted',action='write')
   open(unit=72,file='ddmdq_wf_t2.out',status='unknown',form='formatted',action='write')
   open(unit=73,file='ddmdq_wf_t3.out',status='unknown',form='formatted',action='write')
   open(unit=74,file='ddmdq_ewald.out',status='unknown',form='formatted',action='write')
   open(unit=76,file='ddmdq_elecstic.out',status='unknown',form='formatted',action='write')
 end if

!Gather the different terms in the tensors and print the result
 ABI_ALLOCATE(ddmdq_red,(2,natpert,natpert,nq1grad))

 if (kptopt==3) then

   iq1pert=matom+8
   do iq1grad=1,nq1grad
     iq1dir=q1grad(2,iq1grad)
     do jatpert=1,natpert
       jatom=pert_atdis(1,jatpert)
       jatdir=pert_atdis(2,jatpert)
       do iatpert=1,natpert
         iatom=pert_atdis(1,iatpert)
         iatdir=pert_atdis(2,iatpert)
  
         if (ddmdq_flg(iatpert,jatpert,iq1grad)==1) then

           !Calculate and save the third order energy derivative
           tmpre=ddmdq_qgradhart(re,iatpert,jatpert,iq1grad)+ddmdqwf(re,iatpert,jatpert,iq1grad)+&
         & half*dyewdq(re,iatdir,iatom,jatdir,jatom,iq1dir)
           tmpim=ddmdq_qgradhart(im,iatpert,jatpert,iq1grad)+ddmdqwf(im,iatpert,jatpert,iq1grad)+&
         & half*dyewdq(im,iatdir,iatom,jatdir,jatom,iq1dir)
           d3etot(re,iatdir,iatom,jatdir,jatom,iq1dir,iq1pert)=tmpre
           d3etot(im,iatdir,iatom,jatdir,jatom,iq1dir,iq1pert)=tmpim

           !Calculate and write the q-gradient of the dynamical matrix (twice
           !the Energy derivative, see Gonze and Lee 1997) in the form of first
           !moment of IFC in real space
           ddmdq_red(re,iatpert,jatpert,iq1grad)=-two*tmpim
           ddmdq_red(im,iatpert,jatpert,iq1grad)=two*tmpre

           if (prtvol==1) then
             !Write individual contributions 
             write(71,'(5(i5,4x),2(1x,f20.10))') iatom, iatdir, jatom, jatdir, iq1dir,       &
           & -two*ddmdqwf_t1(:,iatpert,jatpert,iq1grad)

             write(72,'(5(i5,4x),2(1x,f20.10))') iatom, iatdir, jatom, jatdir, iq1dir,       &
           & -two*ddmdqwf_t2(:,iatpert,jatpert,iq1grad)

             write(73,'(5(i5,4x),2(1x,f20.10))') iatom, iatdir, jatom, jatdir, iq1dir,       &
           & -two*ddmdqwf_t3(:,iatpert,jatpert,iq1grad)

             write(74,'(5(i5,4x),2(1x,f20.10))') iatom, iatdir, jatom, jatdir, iq1dir,       &
           & -dyewdq(:,iatdir,iatom,jatdir,jatom,iq1dir)

             write(76,'(5(i5,4x),2(1x,f20.10))') iatom, iatdir, jatom, jatdir, iq1dir,       &
           & -two*ddmdq_qgradhart(:,iatpert,jatpert,iq1grad)
           end if

         end if

       end do
     end do
   end do

 else if (kptopt==2) then

   iq1pert=matom+8
   do iq1grad=1,nq1grad
     iq1dir=q1grad(2,iq1grad)
     do jatpert=1,natpert
       jatom=pert_atdis(1,jatpert)
       jatdir=pert_atdis(2,jatpert)
       do iatpert=1,natpert
         iatom=pert_atdis(1,iatpert)
         iatdir=pert_atdis(2,iatpert)
  
         if (ddmdq_flg(iatpert,jatpert,iq1grad)==1) then

           !Calculate and save the third order energy derivative
           tmpim=ddmdq_qgradhart(im,iatpert,jatpert,iq1grad)+ddmdqwf(im,iatpert,jatpert,iq1grad)+&
         & half*dyewdq(im,iatdir,iatom,jatdir,jatom,iq1dir)
           d3etot(re,iatdir,iatom,jatdir,jatom,iq1dir,iq1pert)=zero
           d3etot(im,iatdir,iatom,jatdir,jatom,iq1dir,iq1pert)=tmpim

           !Calculate and write the q-gradient of the dynamical matrix (twice
           !the Energy derivative, see Gonze and Lee 1997) in the form of first
           !moment of IFC in real space
           ddmdq_red(re,iatpert,jatpert,iq1grad)=-two*tmpim
           ddmdq_red(im,iatpert,jatpert,iq1grad)=zero

           if (prtvol==1) then
             !Write individual contributions 
             write(71,'(5(i5,4x),2(1x,f20.10))') iatom, iatdir, jatom, jatdir, iq1dir,       &
           & -two*ddmdqwf_t1(im,iatpert,jatpert,iq1grad)

             write(72,'(5(i5,4x),2(1x,f20.10))') iatom, iatdir, jatom, jatdir, iq1dir,       &
           & -two*ddmdqwf_t2(im,iatpert,jatpert,iq1grad)

             write(73,'(5(i5,4x),2(1x,f20.10))') iatom, iatdir, jatom, jatdir, iq1dir,       &
           & -two*ddmdqwf_t3(im,iatpert,jatpert,iq1grad)

             write(74,'(5(i5,4x),2(1x,f20.10))') iatom, iatdir, jatom, jatdir, iq1dir,       &
           & -dyewdq(im,iatdir,iatom,jatdir,jatom,iq1dir)

             write(76,'(5(i5,4x),2(1x,f20.10))') iatom, iatdir, jatom, jatdir, iq1dir,       &
           & -two*ddmdq_qgradhart(im,iatpert,jatpert,iq1grad)
           end if

         end if

       end do
     end do
   end do

 else

   write(msg,"(1a)") 'kptopt must be 2 or 3 for the ddmdq calculation'
   MSG_BUG(msg)

 end if

 if (prtvol==1) then
   close(71)
   close(72)
   close(73)
   close(74)
   close(76)
 end if

!Transformation to cartesian coordinates of the ddmdq
 ABI_ALLOCATE(ddmdq_cart,(2,matom,3,matom,3,3))
 ABI_ALLOCATE(ddmdq_cartflg,(matom,3,matom,3,3))
 ABI_ALLOCATE(cartflg,(matom,3,matom,3,3))
 cartflg=0
 do iq1grad=1,nq1grad
   iq1dir=q1grad(2,iq1grad)
   do jatpert=1,natpert
     jatom=pert_atdis(1,jatpert)
     jatdir=pert_atdis(2,jatpert)
     do iatpert=1,natpert
       iatom=pert_atdis(1,iatpert)
       iatdir=pert_atdis(2,iatpert)
       ddmdq_cartflg(iatom,iatdir,jatom,jatdir,iq1dir)=ddmdq_flg(iatpert,jatpert,iq1grad)
       ddmdq_cart(:,iatom,iatdir,jatom,jatdir,iq1dir)=ddmdq_red(:,iatpert,jatpert,iq1grad)
     end do
   end do
 end do
 ABI_DEALLOCATE(ddmdq_red)
 
!1st transform coordenates of the first atomic displacement derivative
 do iq1dir=1,3
   do jatdir=1,3
     do jatom=1,matom
       do ii=1,2
         do iatom=1,matom

           do iatdir=1,3
             vec1(iatdir)=ddmdq_cart(ii,iatom,iatdir,jatom,jatdir,iq1dir)
             flg1(iatdir)=ddmdq_cartflg(iatom,iatdir,jatom,jatdir,iq1dir)
           end do 
           call cart39(flg1,flg2,gprimd,iatom,matom,rprimd,vec1,vec2)
           do iatdir=1,3
             ddmdq_cart(ii,iatom,iatdir,jatom,jatdir,iq1dir)=vec2(iatdir)
             cartflg(iatom,iatdir,jatom,jatdir,iq1dir)=flg2(iatdir)
           end do

         end do
       end do
     end do
   end do
 end do

!2nd transform coordenates of the second atomic displacement derivative
 do iq1dir=1,3
   do iatdir=1,3
     do iatom=1,matom
       do ii=1,2
         do jatom=1,matom

           do jatdir=1,3
             vec1(jatdir)=ddmdq_cart(ii,iatom,iatdir,jatom,jatdir,iq1dir)
             flg1(jatdir)=ddmdq_cartflg(iatom,iatdir,jatom,jatdir,iq1dir)
           end do 
           call cart39(flg1,flg2,gprimd,jatom,matom,rprimd,vec1,vec2)
           do jatdir=1,3
             ddmdq_cart(ii,iatom,iatdir,jatom,jatdir,iq1dir)=vec2(jatdir)
           end do

         end do
       end do
     end do
   end do
 end do

!3rd transform coordinates of the q-gradient (treat it as electric field)
 do jatdir=1,3
   do jatom=1,matom
     do iatdir=1,3
       do iatom=1,matom
         do ii=1,2

           do iq1dir=1,3
             vec1(iq1dir)=ddmdq_cart(ii,iatom,iatdir,jatom,jatdir,iq1dir)
             flg1(iq1dir)=ddmdq_cartflg(iatom,iatdir,jatom,jatdir,iq1dir)
           end do 
           call cart39(flg1,flg2,gprimd,matom+2,matom,rprimd,vec1,vec2)
           do iq1dir=1,3
             ddmdq_cart(ii,iatom,iatdir,jatom,jatdir,iq1dir)=vec2(iq1dir)
           end do

         end do
       end do
     end do
   end do
 end do

!Write the tensor in cartesian coordinates
 write(ab_out,'(a)')' '
 write(ab_out,'(a)')' First moment of real space IFC, in cartesian coordinates,'
 write(ab_out,'(a)')' iatom   iatdir   jatom   jatddir   qgrdir           real part          imaginary part'
 do iq1dir=1,3
   do jatdir=1,3
     do jatom=1,matom
       do iatdir=1,3
         do iatom=1,matom
     
           if (cartflg(iatom,iatdir,jatom,jatdir,iq1dir)==1) then
             write(ab_out,'(5(i5,4x),2(1x,f20.10))') iatom, iatdir, jatom, jatdir, iq1dir,       &
           & ddmdq_cart(:,iatom,iatdir,jatom,jatdir,iq1dir)
           end if


         end do
       end do
     end do
   end do
   write(ab_out,'(a)')' '
 end do
 write(ab_out,'(80a)')('=',ii=1,80)

!Write the piezoelectric force-response tensor
 write(ab_out,'(a)')' '
 write(ab_out,'(a)')' Piezoelectric force-response tensor, in cartesian coordinates '
 write(ab_out,'(a)')' (from q-gradient of dnamical matrix),'
 write(ab_out,'(a)')' iatom   iatddir  jatddir   qgrdir           real part          imaginary part'
 do iq1dir=1,3
   do iatdir=1,3
     do iatom=1,matom
       do jatdir=1,3
         piezofrre=zero
         piezofrim=zero
         do jatom=1,matom

           if (cartflg(iatom,iatdir,jatom,jatdir,iq1dir)==1) then
             piezofrre=piezofrre+ddmdq_cart(1,iatom,iatdir,jatom,jatdir,iq1dir)
             piezofrim=piezofrim+ddmdq_cart(2,iatom,iatdir,jatom,jatdir,iq1dir)
           end if

         end do

         write(ab_out,'(4(i5,4x),2(1x,f20.10))') iatom, iatdir, jatdir, iq1dir,       &
       & piezofrre, piezofrim

       end do
     end do
   end do
   write(ab_out,'(a)')' '
 end do
 write(ab_out,'(80a)')('=',ii=1,80)

 ABI_DEALLOCATE(ddmdq_cart)
 ABI_DEALLOCATE(cartflg)


 DBG_EXIT("COLL")
 end subroutine dfpt_ddmdqout
!!***

!!****f* ABINIT/dfpt_isdqout
!! NAME
!!  dfpt_isdqout
!!
!! FUNCTION
!!  This subroutine gathers the different terms entering the q-gradient of the
!!  internal strain tensor, perfofms the transformation from reduced to cartesian coordinates and 
!!  writes out the tensor in output files in type-II formulation.
!!  
!! COPYRIGHT
!!  Copyright (C) 2018 ABINIT group (MR,MS)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dyewdqdq(2,3,natom,3,3,3)= Second q-gradient of Ewald part of the dynamical matrix
!!           sumed over the second atomic sublattice
!!  frwfdq(2,natpert,3,3,nq1grad)=frozen wf contribution (type-I)
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  isdq_flg(matom,3,nq1grad,3,3)=array that indicates which elements of the q-gradient of
!!                                  internal strain tensor have been calculated
!!  isdq_qgradhart(2,matom,3,nq1grad,3,3)=electronic electrostatic contribution from the 
!!                                             q-gradient of the Hartree potential
!!  isdqwf(2,matom,3,nq1grad,3,3)=total wave function contributions to the q-gradient of
!!                          internal strain tensor (except t4)
!!  isdqwf_t1(2,matom,3,nq1grad,3,3)=term 1 of the wave function contribution 
!!  isdqwf_t2(2,matom,3,nq1grad,3,3)=term 2 of the wave function contribution 
!!  isdqwf_t3(2,matom,3,nq1grad,3,3)=term 3 of the wave function contribution 
!!  isdqwf_t4(2,matom,3,3,3,nq1grad)=term 4 of the wave function contribution (type-I) 
!!  isdqwf_t5(2,matom,3,nq1grad,3,3)=term 5 of the wave function contribution 
!!  kptopt=2 time reversal symmetry is enforced, 3 trs is not enforced (for debugging purposes)
!!  matom=number of atoms 
!!  mpert=maximum number of perturbations
!!  natpert=number of atomic displacement perturbations
!!  nstrpert=number of strain perturbations
!!  nq1grad=number of q1 (q_{\gamma}) gradients
!!  pert_atdis(3,natpert)=array with the info for the atomic displacement perturbations
!!  pert_strain(6,nstrpert)=array with the info for the strain perturbations
!!  prtvol=volume of information to be printed. 1-> The different contributions to the quadrupole are printed.
!!  q1grad(3,nq1grad)=array with the info for the q1 (q_{\gamma}) gradients
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  ucvol=Unit cell volume
!!  
!! OUTPUT
!!  d3etot(2,3,mpert,3,mpert,3,mpert)= matrix of the 3DTE
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!!  dfpt_flexo
!!
!! CHILDREN
!!
!!  cart39 
!!
!! SOURCE

 subroutine dfpt_isdqout(d3etot,dyewdqdq,frwfdq,gprimd,isdq_flg,isdq_qgradhart,isdqwf,isdqwf_t1,isdqwf_t2,&
   & isdqwf_t3,isdqwf_t4,isdqwf_t5,kptopt,matom,mpert,natpert, &
   & nstrpert,nq1grad,pert_atdis,pert_strain,prtvol,q1grad,rprimd,ucvol)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_isdqout'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: kptopt,matom,mpert,natpert,nstrpert,nq1grad,prtvol
 real(dp),intent(in) :: ucvol

!arrays
 integer,intent(in) :: isdq_flg(matom,3,nq1grad,3,3)
 integer,intent(in) :: pert_atdis(3,natpert)
 integer,intent(in) :: pert_strain(6,nstrpert)
 integer,intent(in) :: q1grad(3,nq1grad)
 real(dp),intent(inout) :: d3etot(2,3,mpert,3,mpert,3,mpert)
 real(dp),intent(inout) :: dyewdqdq(2,3,matom,3,3,3)
 real(dp),intent(in) :: isdqwf(2,matom,3,nq1grad,3,3)
 real(dp),intent(inout) :: frwfdq(2,matom,3,3,3,nq1grad)
 real(dp),intent(inout) :: isdq_qgradhart(2,matom,3,nq1grad,3,3)
 real(dp),intent(inout) :: isdqwf_t1(2,matom,3,nq1grad,3,3)
 real(dp),intent(inout) :: isdqwf_t2(2,matom,3,nq1grad,3,3)
 real(dp),intent(inout) :: isdqwf_t3(2,matom,3,nq1grad,3,3)
 real(dp),intent(inout) :: isdqwf_t4(2,matom,3,3,3,nq1grad)
 real(dp),intent(inout) :: isdqwf_t5(2,matom,3,nq1grad,3,3)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: rprimd(3,3)
 
!Local variables-------------------------------
!scalars
 integer :: alpha,beta,delta,gamma
 integer :: iatdir,iatom,iatpert,ibuf,ii,iq1dir,iq1grad,istr1dir,istr2dir,istrpert
 integer :: q1pert,strcomp,strpert
 integer, parameter :: re=1,im=2
 real(dp) :: celastre,celastim,fac,tewim,tewre,tfrim,tfrre,t4im,tmpim,tmpre,t4re
 character(len=500) :: msg                   

!arrays
 integer :: flg1(3),flg2(3)
 integer, allocatable :: typeI_cartflag(:,:,:,:,:),redflg(:,:,:,:,:)
 real(dp) :: vec1(3),vec2(3)
 real(dp),allocatable :: dyewdqdq_cart(:,:,:,:,:,:),frwfdq_cart(:,:,:,:,:,:)
 real(dp),allocatable :: isdqtens_cart(:,:,:,:,:,:),isdqtens_red(:,:,:,:,:,:)
 real(dp),allocatable :: isdqwf_t4_cart(:,:,:,:,:,:)
 real(dp),allocatable :: isdqtens_buffer_cart(:,:,:,:,:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

!Gather the different terms in the q-gradient of the internal strain tensor
 ABI_ALLOCATE(isdqtens_red,(2,matom,3,nq1grad,3,3))
! isdqtens_red=zero

 if (kptopt==3) then

   !Compute real and 'true' imaginary parts of isdq tensor and independent
   !terms. The T4 term and the frozen wf contributions need further treatment 
   !and they will be lately added to the cartesian coordinates
   !version of the isdq tensor
   do istrpert=1,nstrpert
     istr1dir=pert_strain(3,istrpert)
     istr2dir=pert_strain(4,istrpert)
     do iq1grad=1,nq1grad
       iq1dir=q1grad(2,iq1grad)
       do iatpert=1,natpert
         iatom=pert_atdis(1,iatpert)
         iatdir=pert_atdis(2,iatpert)

         if (isdq_flg(iatom,iatdir,iq1dir,istr1dir,istr2dir)==1) then  
        
           !The interesting magnitude is 2 times the gradient of the second order Energy
           isdqtens_red(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)=two* &
         & ( isdq_qgradhart(re,iatom,iatdir,iq1dir,istr1dir,istr2dir) + &
         &   isdqwf(re,iatom,iatdir,iq1dir,istr1dir,istr2dir) )
           isdqtens_red(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)=two* &
         & ( isdq_qgradhart(im,iatom,iatdir,iq1dir,istr1dir,istr2dir) + &
         &   isdqwf(im,iatom,iatdir,iq1dir,istr1dir,istr2dir) )

           !Multiply by the imaginary unit that has been factorized out
           tmpre=isdqtens_red(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)
           tmpim=isdqtens_red(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)
           isdqtens_red(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)=-tmpim
           isdqtens_red(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)=tmpre

           !Do the smae for the T4 term 
           !(This is different from the flexoout routine because the factorized -i in the T4
           ! has to be multiplied by -2 as we did for the rest of contributions)
           tmpre=isdqwf_t4(re,iatom,iatdir,istr1dir,istr2dir,iq1dir)
           tmpim=isdqwf_t4(im,iatom,iatdir,istr1dir,istr2dir,iq1dir)
           isdqwf_t4(re,iatom,iatdir,istr1dir,istr2dir,iq1dir)=two*tmpim
           isdqwf_t4(im,iatom,iatdir,istr1dir,istr2dir,iq1dir)=-two*tmpre

           !Multiply by 2 the frozen wf contribution 
           tmpre=frwfdq(re,iatom,iatdir,istr1dir,istr2dir,iq1dir)
           tmpim=frwfdq(im,iatom,iatdir,istr1dir,istr2dir,iq1dir)
           frwfdq(re,iatom,iatdir,istr1dir,istr2dir,iq1dir)=two*tmpre
           frwfdq(im,iatom,iatdir,istr1dir,istr2dir,iq1dir)=two*tmpim

           !Multiply by 2 the Ewald contribution 
           tmpre=dyewdqdq(re,iatdir,iatom,istr1dir,istr2dir,iq1dir)
           tmpim=dyewdqdq(im,iatdir,iatom,istr1dir,istr2dir,iq1dir)
           dyewdqdq(re,iatdir,iatom,istr1dir,istr2dir,iq1dir)=two*tmpre
           dyewdqdq(im,iatdir,iatom,istr1dir,istr2dir,iq1dir)=two*tmpim
           
           !Compute and save individual terms in mixed coordinates
           if (prtvol==1) then
             
             tmpre=isdq_qgradhart(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)
             tmpim=isdq_qgradhart(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)
             isdq_qgradhart(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)=-two*tmpim
             isdq_qgradhart(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)=two*tmpre

             tmpre=isdqwf_t1(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)
             tmpim=isdqwf_t1(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)
             isdqwf_t1(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)=-two*tmpim
             isdqwf_t1(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)=two*tmpre

             tmpre=isdqwf_t2(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)
             tmpim=isdqwf_t2(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)
             isdqwf_t2(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)=-two*tmpim
             isdqwf_t2(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)=two*tmpre

             tmpre=isdqwf_t3(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)
             tmpim=isdqwf_t3(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)
             isdqwf_t3(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)=-two*tmpim
             isdqwf_t3(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)=two*tmpre

             tmpre=isdqwf_t5(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)
             tmpim=isdqwf_t5(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)
             isdqwf_t5(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)=-two*tmpim
             isdqwf_t5(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)=two*tmpre

           end if

         end if

       end do
     end do
   end do

 else if (kptopt==2) then

   !Compute real part of isdq tensor and independent
   !terms. The T4 term and the frozen wf contributions need further treatment 
   !and they will be lately added to the cartesian coordinates
   !version of the isdq tensor
   do istrpert=1,nstrpert
     istr1dir=pert_strain(3,istrpert)
     istr2dir=pert_strain(4,istrpert)
     do iq1grad=1,nq1grad
       iq1dir=q1grad(2,iq1grad)
       do iatpert=1,natpert
         iatom=pert_atdis(1,iatpert)
         iatdir=pert_atdis(2,iatpert)

         if (isdq_flg(iatom,iatdir,iq1dir,istr1dir,istr2dir)==1) then  
        
           !The interesting magnitude is 2 times the gradient of the second order Energy
           isdqtens_red(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)=two* &
         & ( isdq_qgradhart(re,iatom,iatdir,iq1dir,istr1dir,istr2dir) + &
         &   isdqwf(re,iatom,iatdir,iq1dir,istr1dir,istr2dir) )
           isdqtens_red(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)=two* &
         & ( isdq_qgradhart(im,iatom,iatdir,iq1dir,istr1dir,istr2dir) + &
         &   isdqwf(im,iatom,iatdir,iq1dir,istr1dir,istr2dir) )

           !Multiply by the imaginary unit that has been factorized out
           tmpre=isdqtens_red(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)
           tmpim=isdqtens_red(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)
           isdqtens_red(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)=-tmpim
           isdqtens_red(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)=zero

           !Do the smae for the T4 term 
           !(This is different from the flexoout routine because the factorized -i in the T4
           ! has to be multiplied by 2 as we did for the rest of contributions)
           tmpre=isdqwf_t4(re,iatom,iatdir,istr1dir,istr2dir,iq1dir)
           tmpim=isdqwf_t4(im,iatom,iatdir,istr1dir,istr2dir,iq1dir)
           isdqwf_t4(re,iatom,iatdir,istr1dir,istr2dir,iq1dir)=two*tmpim
           isdqwf_t4(im,iatom,iatdir,istr1dir,istr2dir,iq1dir)=zero

           !Multiply by 2 the frozen wf contribution 
           tmpre=frwfdq(re,iatom,iatdir,istr1dir,istr2dir,iq1dir)
           frwfdq(re,iatom,iatdir,istr1dir,istr2dir,iq1dir)=two*tmpre
           frwfdq(im,iatom,iatdir,istr1dir,istr2dir,iq1dir)=zero

           !Multiply by 2 the Ewald contribution 
           tmpre=dyewdqdq(re,iatdir,iatom,istr1dir,istr2dir,iq1dir)
           dyewdqdq(re,iatdir,iatom,istr1dir,istr2dir,iq1dir)=two*tmpre
           dyewdqdq(im,iatdir,iatom,istr1dir,istr2dir,iq1dir)=zero
           
           !Compute and save individual terms in mixed coordinates
           if (prtvol==1) then
             
             tmpim=isdq_qgradhart(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)
             isdq_qgradhart(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)=-two*tmpim
             isdq_qgradhart(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)=zero

             tmpim=isdqwf_t1(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)
             isdqwf_t1(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)=-two*tmpim
             isdqwf_t1(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)=zero

             tmpim=isdqwf_t2(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)
             isdqwf_t2(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)=-two*tmpim
             isdqwf_t2(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)=zero

             tmpim=isdqwf_t3(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)
             isdqwf_t3(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)=-two*tmpim
             isdqwf_t3(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)=zero

             tmpim=isdqwf_t5(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)
             isdqwf_t5(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)=-two*tmpim
             isdqwf_t5(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)=zero

           end if

         end if

       end do
     end do
   end do

 else
   write(msg,"(1a)") 'kptopt must be 2 or 3 for long-wave DFPT calculations'
   MSG_BUG(msg)
 end if

!Transform to complete cartesian coordinates all the contributions
 ABI_ALLOCATE(isdqtens_cart,(2,matom,3,nq1grad,3,3))
 ABI_ALLOCATE(isdqwf_t4_cart,(2,matom,3,3,3,nq1grad))
 ABI_ALLOCATE(frwfdq_cart,(2,matom,3,3,3,nq1grad))
 ABI_ALLOCATE(dyewdqdq_cart,(2,3,matom,3,3,nq1grad))
 ABI_ALLOCATE(typeI_cartflag,(matom,3,3,3,nq1grad))
 isdqtens_cart=isdqtens_red
 isdqwf_t4_cart=isdqwf_t4
 frwfdq_cart=frwfdq
 dyewdqdq_cart=dyewdqdq
 typeI_cartflag=0

 if (prtvol==1) then
   ABI_ALLOCATE(isdqtens_buffer_cart,(5,2,matom,3,nq1grad,3,3))
   isdqtens_buffer_cart(1,:,:,:,:,:,:)=isdqwf_t1(:,:,:,:,:,:)
   isdqtens_buffer_cart(2,:,:,:,:,:,:)=isdqwf_t2(:,:,:,:,:,:)
   isdqtens_buffer_cart(3,:,:,:,:,:,:)=isdqwf_t3(:,:,:,:,:,:)
   isdqtens_buffer_cart(4,:,:,:,:,:,:)=isdq_qgradhart(:,:,:,:,:,:)
   isdqtens_buffer_cart(5,:,:,:,:,:,:)=isdqwf_t5(:,:,:,:,:,:)
 end if

!##1st transform coordinates of the atomic displacement derivative contributions
 do istr2dir=1,3
   do istr1dir=1,3
     do iq1dir=1,3
       do ii=1,2
         do iatom=1,matom

           do iatdir=1,3
             vec1(iatdir)=isdqtens_cart(ii,iatom,iatdir,iq1dir,istr1dir,istr2dir)
             flg1(iatdir)=isdq_flg(iatom,iatdir,iq1dir,istr1dir,istr2dir)
           end do
           call cart39(flg1,flg2,gprimd,iatom,matom,rprimd,vec1,vec2)
           do iatdir=1,3
             isdqtens_cart(ii,iatom,iatdir,iq1dir,istr1dir,istr2dir)=vec2(iatdir)
           end do

           do iatdir=1,3
             vec1(iatdir)=isdqwf_t4_cart(ii,iatom,iatdir,istr1dir,istr2dir,iq1dir)
             flg1(iatdir)=isdq_flg(iatom,iatdir,iq1dir,istr1dir,istr2dir)
           end do
           call cart39(flg1,flg2,gprimd,iatom,matom,rprimd,vec1,vec2)
           do iatdir=1,3
             isdqwf_t4_cart(ii,iatom,iatdir,istr1dir,istr2dir,iq1dir)=vec2(iatdir)
             typeI_cartflag(iatom,iatdir,istr1dir,istr2dir,iq1dir)=flg2(iatdir)
           end do

           do iatdir=1,3
             vec1(iatdir)=frwfdq_cart(ii,iatom,iatdir,istr1dir,istr2dir,iq1dir)
             flg1(iatdir)=isdq_flg(iatom,iatdir,iq1dir,istr1dir,istr2dir)
           end do
           call cart39(flg1,flg2,gprimd,iatom,matom,rprimd,vec1,vec2)
           do iatdir=1,3
             frwfdq_cart(ii,iatom,iatdir,istr1dir,istr2dir,iq1dir)=vec2(iatdir)
           end do

           do iatdir=1,3
             vec1(iatdir)=dyewdqdq_cart(ii,iatdir,iatom,istr1dir,istr2dir,iq1dir)
             flg1(iatdir)=isdq_flg(iatom,iatdir,iq1dir,istr1dir,istr2dir)
           end do
           call cart39(flg1,flg2,gprimd,iatom,matom,rprimd,vec1,vec2)
           do iatdir=1,3
             dyewdqdq_cart(ii,iatdir,iatom,istr1dir,istr2dir,iq1dir)=vec2(iatdir)
           end do

         end do 
       end do
     end do
   end do
 end do

!For debugging, transform also all the individual contributions
 if (prtvol==1) then
   do ibuf=1,5
     do istr2dir=1,3
       do istr1dir=1,3
         do iq1dir=1,3
           do ii=1,2
             do iatom=1,matom
    
               do iatdir=1,3
                 vec1(iatdir)=isdqtens_buffer_cart(ibuf,ii,iatom,iatdir,iq1dir,istr1dir,istr2dir)
                 flg1(iatdir)=isdq_flg(iatom,iatdir,iq1dir,istr1dir,istr2dir)
               end do
               call cart39(flg1,flg2,gprimd,iatom,matom,rprimd,vec1,vec2)
               do iatdir=1,3
                 isdqtens_buffer_cart(ibuf,ii,iatom,iatdir,iq1dir,istr1dir,istr2dir)=vec2(iatdir)
               end do

             end do 
           end do
         end do
       end do
     end do
   end do
 end if

!##2nd transform coordinates of the q-gradient (treat it as electric field)
 do istr2dir=1,3
   do istr1dir=1,3
     do iatom=1,matom
       do iatdir=1,3
         do ii=1,2

           do iq1dir=1,3
             vec1(iq1dir)=isdqtens_cart(ii,iatom,iatdir,iq1dir,istr1dir,istr2dir)
             flg1(iq1dir)=isdq_flg(iatom,iatdir,iq1dir,istr1dir,istr2dir)
           end do
           call cart39(flg1,flg2,gprimd,matom+2,matom,rprimd,vec1,vec2)
           do iq1dir=1,3
             isdqtens_cart(ii,iatom,iatdir,iq1dir,istr1dir,istr2dir)=vec2(iq1dir)
           end do

           do iq1dir=1,3
             vec1(iq1dir)=frwfdq_cart(ii,iatom,iatdir,istr1dir,istr2dir,iq1dir)
             flg1(iq1dir)=isdq_flg(iatom,iatdir,iq1dir,istr1dir,istr2dir)
           end do
           call cart39(flg1,flg2,gprimd,matom+2,matom,rprimd,vec1,vec2)
           do iq1dir=1,3
             frwfdq_cart(ii,iatom,iatdir,istr1dir,istr2dir,iq1dir)=vec2(iq1dir)
           end do

           do iq1dir=1,3
             vec1(iq1dir)=dyewdqdq_cart(ii,iatdir,iatom,istr1dir,istr2dir,iq1dir)
             flg1(iq1dir)=isdq_flg(iatom,iatdir,iq1dir,istr1dir,istr2dir)
           end do
           call cart39(flg1,flg2,gprimd,matom+2,matom,rprimd,vec1,vec2)
           do iq1dir=1,3
             dyewdqdq_cart(ii,iatdir,iatom,istr1dir,istr2dir,iq1dir)=vec2(iq1dir)
           end do

         end do 
       end do
     end do
   end do
 end do

!For debugging, transform also all the individual contributions
 if (prtvol==1) then
   do ibuf=1,5
     do istr2dir=1,3
       do istr1dir=1,3
         do iatom=1,matom
           do iatdir=1,3
             do ii=1,2
    
               do iq1dir=1,3
                 vec1(iq1dir)=isdqtens_buffer_cart(ibuf,ii,iatom,iatdir,iq1dir,istr1dir,istr2dir)
                 flg1(iq1dir)=isdq_flg(iatom,iatdir,iq1dir,istr1dir,istr2dir)
               end do
               call cart39(flg1,flg2,gprimd,matom+2,matom,rprimd,vec1,vec2)
               do iq1dir=1,3
                 isdqtens_buffer_cart(ibuf,ii,iatom,iatdir,iq1dir,istr1dir,istr2dir)=vec2(iq1dir)
               end do
    
             end do 
           end do
         end do
       end do
     end do
   end do
 endif

!3rd## Transform the metric perturbation direction to cartesian coordinates in
!the frozen wf contributions (treat it as an atomic displacement)
 do istr2dir=1,3
   do iq1dir=1,3
     do iatom=1,matom
       do iatdir=1,3
         do ii=1,2

           do istr1dir=1,3
             vec1(istr1dir)=frwfdq_cart(ii,iatom,iatdir,istr1dir,istr2dir,iq1dir)
             flg1(istr1dir)=isdq_flg(iatom,iatdir,iq1dir,istr1dir,istr2dir)
           end do
           call cart39(flg1,flg2,gprimd,iatom,matom,rprimd,vec1,vec2)
           do istr1dir=1,3
             frwfdq_cart(ii,iatom,iatdir,istr1dir,istr2dir,iq1dir)=vec2(istr1dir)
           end do

           do istr1dir=1,3
             vec1(istr1dir)=dyewdqdq_cart(ii,iatdir,iatom,istr1dir,istr2dir,iq1dir)
             flg1(istr1dir)=isdq_flg(iatom,iatdir,iq1dir,istr1dir,istr2dir)
           end do
           call cart39(flg1,flg2,gprimd,iatom,matom,rprimd,vec1,vec2)
           do istr1dir=1,3
             dyewdqdq_cart(ii,iatdir,iatom,istr1dir,istr2dir,iq1dir)=vec2(istr1dir)
           end do

         end do 
       end do
     end do
   end do
 end do

!4th## Transform the second q-gradient direction to cartesian coordinates in
!the frozen wf contributions (treat it as an electric field)
 do istr1dir=1,3
   do iq1dir=1,3
     do iatom=1,matom
       do iatdir=1,3
         do ii=1,2

           do istr2dir=1,3
             vec1(istr2dir)=frwfdq_cart(ii,iatom,iatdir,istr1dir,istr2dir,iq1dir)
             flg1(istr2dir)=isdq_flg(iatom,iatdir,iq1dir,istr1dir,istr2dir)
           end do
           call cart39(flg1,flg2,gprimd,matom+2,matom,rprimd,vec1,vec2)
           do istr2dir=1,3
             frwfdq_cart(ii,iatom,iatdir,istr1dir,istr2dir,iq1dir)=vec2(istr2dir)
           end do

           do istr2dir=1,3
             vec1(istr2dir)=dyewdqdq_cart(ii,iatdir,iatom,istr1dir,istr2dir,iq1dir)
             flg1(istr2dir)=isdq_flg(iatom,iatdir,iq1dir,istr1dir,istr2dir)
           end do
           call cart39(flg1,flg2,gprimd,matom+2,matom,rprimd,vec1,vec2)
           do istr2dir=1,3
             dyewdqdq_cart(ii,iatdir,iatom,istr1dir,istr2dir,iq1dir)=vec2(istr2dir)
           end do

         end do 
       end do
     end do
   end do
 end do

!Write the q-gradient fo the internal strain tensor in cartesian coordinates
!Open output files
 if (prtvol==1) then
   open(unit=71,file='isdq_wf_t1.out',status='unknown',form='formatted',action='write')
   open(unit=72,file='isdq_wf_t2.out',status='unknown',form='formatted',action='write')
   open(unit=73,file='isdq_wf_t3.out',status='unknown',form='formatted',action='write')
   open(unit=74,file='isdq_wf_t4.out',status='unknown',form='formatted',action='write')
   open(unit=75,file='isdq_wf_t5.out',status='unknown',form='formatted',action='write')
   open(unit=76,file='isdq_elecstic.out',status='unknown',form='formatted',action='write')
   open(unit=77,file='isdq_frwf.out',status='unknown',form='formatted',action='write')
   open(unit=78,file='isdq_ewdqdq.out',status='unknown',form='formatted',action='write')
 end if

 write(ab_out,'(a)')' '
 write(ab_out,'(a)')' q-gradient of piezoelectric force response tensor, in cartesian coordinates,'
 write(ab_out,'(a)')'atom   atdir  qgrdir  strdir1  strdir2         real part          imaginary part'
 do istr2dir=1,3
   delta=istr2dir
   do istr1dir=1,3
     beta=istr1dir
     do iq1dir=1,3
       gamma=iq1dir
       do iatom=1,matom
         do iatdir=1,3
           alpha=iatdir

           if (typeI_cartflag(iatom,alpha,beta,delta,gamma)==1 .and. &
          &    typeI_cartflag(iatom,alpha,delta,gamma,beta)==1 .and. & 
          &    typeI_cartflag(iatom,alpha,gamma,beta,delta)==1 ) then

             !Converts the T4 term to type-II form
             t4re= isdqwf_t4_cart(re,iatom,alpha,beta,delta,gamma) + &
                 & isdqwf_t4_cart(re,iatom,alpha,delta,gamma,beta) - &
                 & isdqwf_t4_cart(re,iatom,alpha,gamma,beta,delta) 
             t4im= isdqwf_t4_cart(im,iatom,alpha,beta,delta,gamma) + &
                 & isdqwf_t4_cart(im,iatom,alpha,delta,gamma,beta) - &
                 & isdqwf_t4_cart(im,iatom,alpha,gamma,beta,delta) 

             !Converts the frozen wf term to type-II form
             tfrre= frwfdq_cart(re,iatom,alpha,beta,delta,gamma) + &
                  & frwfdq_cart(re,iatom,alpha,delta,gamma,beta) - & 
                  & frwfdq_cart(re,iatom,alpha,gamma,beta,delta) 
             tfrim= frwfdq_cart(im,iatom,alpha,beta,delta,gamma) + &
                  & frwfdq_cart(im,iatom,alpha,delta,gamma,beta) - & 
                  & frwfdq_cart(im,iatom,alpha,gamma,beta,delta) 

             !Converts the Ewald term to type-II form
             tewre= dyewdqdq_cart(re,alpha,iatom,beta,delta,gamma) + &
                  & dyewdqdq_cart(re,alpha,iatom,delta,gamma,beta) - & 
                  & dyewdqdq_cart(re,alpha,iatom,gamma,beta,delta) 

             tewim= dyewdqdq_cart(im,alpha,iatom,beta,delta,gamma) + &
                  & dyewdqdq_cart(im,alpha,iatom,delta,gamma,beta) - & 
                  & dyewdqdq_cart(im,alpha,iatom,gamma,beta,delta) 

             !Add the type-II T4, frozen wf and Ewald contributions 
             isdqtens_cart(re,iatom,alpha,gamma,beta,delta)= &
           & isdqtens_cart(re,iatom,alpha,gamma,beta,delta) + t4re + &
           & half*tfrre + quarter*tewre
             isdqtens_cart(im,iatom,alpha,gamma,beta,delta)= &
           & isdqtens_cart(im,iatom,alpha,gamma,beta,delta) + t4im + &
           & half*tfrim + quarter*tewim

             !Writes the complete q-gradient of internal strain tensor
             write(ab_out,'(5(i5,3x),2(1x,f20.10))') iatom,alpha,gamma,beta,delta, &
           & isdqtens_cart(re,iatom,alpha,gamma,beta,delta), &
           & isdqtens_cart(im,iatom,alpha,gamma,beta,delta)
 
             if (prtvol==1) then
               write(71,'(5(i5,3x),2(1x,f20.10))') iatom,alpha,gamma,beta,delta, &
             & isdqtens_buffer_cart(1,re,iatom,alpha,gamma,beta,delta), &
             & isdqtens_buffer_cart(1,im,iatom,alpha,gamma,beta,delta)
   
               write(72,'(5(i5,3x),2(1x,f20.10))') iatom,alpha,gamma,beta,delta, &
             & isdqtens_buffer_cart(2,re,iatom,alpha,gamma,beta,delta), &
             & isdqtens_buffer_cart(2,im,iatom,alpha,gamma,beta,delta)

               write(73,'(5(i5,3x),2(1x,f20.10))') iatom,alpha,gamma,beta,delta, &
             & isdqtens_buffer_cart(3,re,iatom,alpha,gamma,beta,delta), &
             & isdqtens_buffer_cart(3,im,iatom,alpha,gamma,beta,delta)

               write(74,'(5(i5,3x),2(1x,f20.10))') iatom,alpha,gamma,beta,delta,t4re,t4im

               write(75,'(5(i5,3x),2(1x,f20.10))') iatom,alpha,gamma,beta,delta, &
             & isdqtens_buffer_cart(5,re,iatom,alpha,gamma,beta,delta), &
             & isdqtens_buffer_cart(5,im,iatom,alpha,gamma,beta,delta)

               write(76,'(5(i5,3x),2(1x,f20.10))') iatom,alpha,gamma,beta,delta, &
             & isdqtens_buffer_cart(4,re,iatom,alpha,gamma,beta,delta), &
             & isdqtens_buffer_cart(4,im,iatom,alpha,gamma,beta,delta)

               write(77,'(5(i5,3x),2(1x,f20.10))') iatom,alpha,gamma,beta,delta, &
             & half*tfrre, half*tfrim

               write(78,'(5(i5,3x),2(1x,f20.10))') iatom,alpha,gamma,beta,delta, &
             & quarter*tewre, quarter*tewim

             end if 

           end if

         end do
       end do
       write(ab_out,'(a)')' '
       if (prtvol==1) then
         write(71,'(a)')' ' 
         write(72,'(a)')' ' 
         write(73,'(a)')' ' 
         write(74,'(a)')' ' 
         write(75,'(a)')' ' 
         write(76,'(a)')' ' 
         write(77,'(a)')' ' 
         write(78,'(a)')' ' 
       end if
     end do
   end do
 end do

 if (prtvol==1) then
   close(71)
   close(72)
   close(73)
   close(74)
   close(75)
   close(76)
   close(77)
   close(78)
 end if

!Write the elastic tensor
 write(ab_out,'(a)')' '
 write(ab_out,'(a)')' Elastic tensor, in cartesian coordinates '
 write(ab_out,'(a)')' (from q-gradient of internal strain),'
 write(ab_out,'(a)')' (for stressed cells it lacks an improper contribution),'
 write(ab_out,'(a)')' atdir  qgrdir  strdir1  strdir2         real part          imaginary part'
 do istr2dir=1,3
   delta=istr2dir
   do istr1dir=1,3
     beta=istr1dir
     do iq1dir=1,3
       gamma=iq1dir
       do iatdir=1,3
         alpha=iatdir
         celastre=zero
         celastim=zero
         do iatom=1,matom
           
           if (isdq_flg(iatom,alpha,gamma,beta,delta)==1) then
             celastre=celastre+isdqtens_cart(re,iatom,alpha,gamma,beta,delta)
             celastim=celastim+isdqtens_cart(im,iatom,alpha,gamma,beta,delta)
           end if

         end do

         write(ab_out,'(4(i5,4x),2(1x,f20.10))') alpha,gamma,beta,delta, &
       & celastre/ucvol, celastim/ucvol

       end do
     end do
   write(ab_out,'(a)')' '
   end do
 end do
 write(ab_out,'(80a)')('=',ii=1,80)

!Calculate the contribution to the d3etot in mixed (reduced/cartesian) coordinates
 isdqtens_red=isdqtens_cart
 ABI_DEALLOCATE(isdqtens_cart) 
 ABI_DEALLOCATE(isdqwf_t4_cart)
 ABI_DEALLOCATE(frwfdq_cart)
 ABI_DEALLOCATE(dyewdqdq_cart)
 ABI_DEALLOCATE(typeI_cartflag)
 ABI_ALLOCATE(redflg,(matom,3,3,3,3))

!1st transform back coordinates of the atomic displacement derivative 
 do istr2dir=1,3
   do istr1dir=1,3
     do iq1dir=1,3
       do ii=1,2
         do iatom=1,matom
           do iatdir=1,3
             vec1(iatdir)=isdqtens_red(ii,iatom,iatdir,iq1dir,istr1dir,istr2dir)
             flg1(iatdir)=isdq_flg(iatom,iatdir,iq1dir,istr1dir,istr2dir)
           end do
           call cart39(flg1,flg2,transpose(rprimd),iatom,matom,transpose(gprimd),vec1,vec2)
           do iatdir=1,3
             isdqtens_red(ii,iatom,iatdir,iq1dir,istr1dir,istr2dir)=vec2(iatdir)
             redflg(iatom,iatdir,iq1dir,istr1dir,istr2dir)=flg2(iatdir)
           end do
         end do
       end do
     end do
   end do
 end do

!2nd transform back coordinates of the q-gradient (treat it as electric field)
 fac=two_pi ** 2
 do istr2dir=1,3
   do istr1dir=1,3
     do iatom=1,matom
       do iatdir=1,3
         do ii=1,2
           do iq1dir=1,3
           vec1(iq1dir)=isdqtens_red(ii,iatom,iatdir,iq1dir,istr1dir,istr2dir)
           flg1(iq1dir)=isdq_flg(iatom,iatdir,iq1dir,istr1dir,istr2dir)
           end do
           call cart39(flg1,flg2,transpose(rprimd),matom+2,matom,transpose(gprimd),vec1,vec2)
           do iq1dir=1,3
             isdqtens_red(ii,iatom,iatdir,iq1dir,istr1dir,istr2dir)=vec2(iq1dir)*fac
             redflg(iatom,iatdir,iq1dir,istr1dir,istr2dir)=flg2(iq1dir)
           end do
         end do
       end do
     end do
   end do
 end do

!Add contributions to d3etot (remove the -2 factor)
 q1pert=matom+8
 do istrpert=1,nstrpert
   strpert=pert_strain(1,istrpert)
   strcomp=pert_strain(2,istrpert)
   istr1dir=pert_strain(3,istrpert)
   istr2dir=pert_strain(4,istrpert)
   do iq1grad=1,nq1grad
     iq1dir=q1grad(2,iq1grad)
     do iatom=1,matom
       do iatdir=1,3

         if (redflg(iatom,iatdir,iq1dir,istr1dir,istr2dir)==1) then
           d3etot(re,iatdir,iatom,strcomp,strpert,iq1dir,q1pert)= &
         & -half*isdqtens_red(re,iatom,iatdir,iq1dir,istr1dir,istr2dir)
           d3etot(im,iatdir,iatom,strcomp,strpert,iq1dir,q1pert)= &
         & -half*isdqtens_red(im,iatom,iatdir,iq1dir,istr1dir,istr2dir)
         end if

       end do
     end do
   end do
 end do

 ABI_DEALLOCATE(isdqtens_red)
 ABI_DEALLOCATE(redflg)

 DBG_EXIT("COLL")
 end subroutine dfpt_isdqout
!!***

end module m_dfpt_lw
!!***
