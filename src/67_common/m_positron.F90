!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_positron
!! NAME
!!  m_positron
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (GJ, MT, JW)
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

module m_positron

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_efield
 use m_errors
 use m_abicore
 use m_energies
 use m_wffile
 use m_electronpositron
 use m_hdr
 use m_xmpi
 use m_bandfft_kpt

 use m_special_funcs,  only : sbf8
 use m_ioarr,    only : ioarr, read_rhor
 use m_pawang,   only : pawang_type, realgaunt
 use m_pawrad,   only : pawrad_type, simp_gen, nderiv_gen
 use m_pawtab,   only : pawtab_type
 use m_paw_ij,   only : paw_ij_type
 use m_pawfgrtab,only : pawfgrtab_type
 use m_pawrhoij,only : pawrhoij_type, pawrhoij_copy, pawrhoij_alloc, pawrhoij_free,&
                       pawrhoij_nullify, pawrhoij_gather, pawrhoij_inquire_dim, pawrhoij_symrhoij
 use m_pawcprj,  only : pawcprj_type, pawcprj_alloc, pawcprj_get, pawcprj_mpi_send, &
                        pawcprj_mpi_recv, pawcprj_free, pawcprj_copy, pawcprj_bcast
 use m_pawfgr,   only : pawfgr_type
 use m_paw_nhat, only : pawmknhat
 use m_fock,     only : fock_type
 use m_kg,       only : getcut
 use defs_wvltypes,     only : wvl_data
 use m_spacepar,        only : hartre
 use m_mkrho,           only : initro
 use m_paw_occupancies, only : initrhoij, pawaccrhoij
 use m_gammapositron, only : gammapositron, gammapositron_fft
 use m_forstr,          only : forstr
 use m_pawxc,         only : pawxcsum
 use m_paw_denpot,    only : pawdensities
 use m_drivexc,       only : mkdenpos

 use m_paw_sphharm, only : initylmr
 use m_pawpsp,  only : pawpsp_read_corewf
 use m_crystal, only : crystal_t
 use m_mpinfo,  only : ptabs_fourdp,set_mpi_enreg_fft,unset_mpi_enreg_fft,destroy_mpi_enreg, initmpi_seq, proc_distrb_cycle
 use m_io_tools,only : open_file,close_unit,get_unit
 use m_fftcore, only : sphereboundary
 use m_prep_kgb,        only : prep_fourwf
 use m_fft,            only : fourwf, fourdp

 implicit none

 private
!!***

 public :: setup_positron
 public :: poslifetime
 public :: posdoppler
!!***

contains
!!***

!!****f* ABINIT/setup_positron
!! NAME
!! setup_positron
!!
!! FUNCTION
!! Do various initializations for the positron lifetime calculation
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  ecore=core psp energy (part of total energy) (hartree)
!!  etotal=current value of total energy
!!  fock <type(fock_type)>= quantities to calculate Fock exact exchange
!!  forces_needed=if >0 forces are needed
!!  fred(3,natom)=forces in reduced coordinates
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space
!!  gmet(3,3)=reciprocal space metric
!!  grchempottn(3,natom)=d(E_chemical_potential)/d(xred) (hartree)
!!  grewtn(3,natom)=d(Ewald)/d(xred) (hartree)
!!  grvdw(3,ngrvdw)=gradients of energy due to Van der Waals DFT-D dispersion (hartree)
!!  gsqcut=cutoff value on G**2 for sphere inside fft box
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  ifirst_gs= 0 if we are in a single ground-state calculation
!!     or in the first ground-state calculation of a structural minimization/dynamics
!!  indsym(4,nsym,natom)=index showing transformation of atom labels
!!                       under symmetry operations (computed in symatm)
!!  istep=index of the number of steps in the routine scfcv
!!  istep_mix=index of the number of steps for the SCF mixing (can be <istep)
!!  kg(3,mpw*mkmem)=reduced (integer) coordinates of G vecs in basis sphere
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only if nkxc>0
!!  maxfor=maximum absolute value of fcart (forces)
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  n3xccc=dimension of the xccc3d array (0 or nfftf).
!!  nattyp(ntypat)= # atoms of each type.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT
!!  ngrvdw=size of grvdw(:,:); can be 0 or natom according to dtset%vdw_xc
!!  nhat(nfftf,nspden*usepaw)= -PAW only- compensation density
!!  nkxc=second dimension of the array kxc, see rhotoxc.f for a description
!!  npwarr(nkpt)=number of planewaves in basis and on boundary for each k
!!  nvresid(nfftf,nspden)=array for the residual of the density/potential
!!  optres=0 if the potential residual has to be used for forces corrections
!!        =1 if the density residual has to be used for forces corrections
!!  paw_ij(my_natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr(my_natom*usepaw) <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(my_natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfftf+1)*natom)=one-dimensional structure factor information (fine FFT grid)
!!  ph1dc(2,3*(2*mgfft+1)*natom)=1-dim structure factor phases (coarse FFT grid)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  stress_needed=if >0 stresses are needed
!!  strsxc(6)=xc correction to stress
!!  symrec(3,3,nsym)=symmetries in reciprocal space, reduced coordinates
!!  ucvol=unit cell volume in bohr**3.
!!  usecprj= 1 if cprj array is stored in memory
!!  vhartr(nfftf)=array for holding Hartree potential
!!  vpsp(nfftf)=array for holding local psp
!!  vxc(nfftf,nspden)=exchange-correlation potential (hartree) in real space
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  energies <type(energies_type)>=all part of total energy.
!!  cg(2,mcg)=wavefunctions
!!  cprj(natom,mcprj*usecprj)= wave functions projected with non-local projectors:
!!                             cprj(n,k,i)=<p_i|Cnk> where p_i is a non-local projector.
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  occ(mband*nkpt*nsppol)=occupation number for each band at each k point
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  rhog(2,nfft)=Fourier transform of total electron/positron density
!!  rhor(nfft,nspden)=total electron/positron density (el/bohr**3)
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      energies_copy,energies_init,forstr,fourdp,getcut,hartre,initrhoij
!!      initro,pawcprj_alloc,pawcprj_copy,pawcprj_free,pawmknhat,pawrhoij_alloc
!!      pawrhoij_copy,pawrhoij_free,read_rhor,wrtout
!!
!! SOURCE

subroutine setup_positron(atindx,atindx1,cg,cprj,dtefield,dtfil,dtset,ecore,eigen,etotal,electronpositron,&
&          energies,fock,forces_needed,fred,gmet,gprimd,grchempottn,grewtn,grvdw,gsqcut,hdr,ifirst_gs,indsym,istep,istep_mix,kg,&
&          kxc,maxfor,mcg,mcprj,mgfft,mpi_enreg,my_natom,n3xccc,nattyp,nfft,ngfft,ngrvdw,nhat,nkxc,npwarr,nvresid,occ,optres,&
&          paw_ij,pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,ph1d,ph1dc,psps,rhog,rhor,&
&          rprimd,stress_needed,strsxc,symrec,ucvol,usecprj,vhartr,vpsp,vxc,&
&          xccc3d,xred,ylm,ylmgr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: forces_needed,ifirst_gs,istep,mcg,mcprj,mgfft,my_natom,n3xccc,nfft
 integer,intent(in) :: ngrvdw,nkxc,optres,stress_needed,usecprj
 integer,intent(inout) :: istep_mix
 real(dp),intent(in) :: ecore,etotal,gsqcut,maxfor,ucvol
 type(efield_type),intent(in) :: dtefield
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer :: electronpositron
 type(energies_type),intent(inout) :: energies
 type(hdr_type),intent(inout) :: hdr
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type), intent(in) :: psps
type(fock_type),pointer, intent(inout) :: fock
!arrays
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom),indsym(4,dtset%nsym,dtset%natom)
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),nattyp(dtset%natom),ngfft(18)
 integer,intent(in) :: npwarr(dtset%nkpt),symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),grchempottn(3,dtset%natom)
 real(dp),intent(in) :: grewtn(3,dtset%natom),grvdw(3,ngrvdw),kxc(nfft,nkxc)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom),ph1dc(2,(3*(2*dtset%mgfft+1)*dtset%natom)*dtset%usepaw)
 real(dp),intent(in) :: rprimd(3,3),strsxc(6),vhartr(nfft),vpsp(nfft),vxc(nfft,dtset%nspden)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(inout) :: cg(2,mcg)
 real(dp),intent(inout) :: nhat(nfft,dtset%nspden*dtset%usepaw)
 real(dp),intent(inout) :: nvresid(nfft,dtset%nspden)
 real(dp),intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol),fred(3,dtset%natom)
 real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(inout) :: rhog(2,nfft),rhor(nfft,dtset%nspden)
 real(dp),intent(inout) :: xccc3d(n3xccc),xred(3,dtset%natom)
 type(pawcprj_type),intent(inout) :: cprj(dtset%natom,mcprj*usecprj)
 type(paw_ij_type),intent(in) :: paw_ij(my_natom*dtset%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*dtset%usepaw)
 type(pawrad_type),intent(in)  :: pawrad(dtset%ntypat*dtset%usepaw)
 type(pawtab_type),intent(in)  :: pawtab(dtset%ntypat*dtset%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: cplex1=1
 integer :: history_level,iatom,iband,icalctype,icalctype0,icg,ifft,ikpt
 integer :: iocc,ireadwf,ispden,isppol,n3xccc0,nocc,optfor,optstr,rdwrpaw,comm_cell
 logical,parameter :: always_restart=.false.  ! Set to true to restart by a pure electronic step at each new atomic structure
 logical :: need_scocc,new_calctype
 real(dp) :: boxcut_dum,diffor_dum,ecut_eff,eigtmp,etotal_read,gsqcut_eff,maxfor_dum
 real(dp) :: maxocc,nelect,occlast,occtmp,rhotmp
 character(len=69) :: TypeCalcStrg
 character(len=500) :: message
 character(len=fnlen) :: fname
 type(energies_type) :: energies_tmp
 type(wvl_data) :: wvl
 type(hdr_type) :: hdr_den
!arrays
 integer,allocatable :: nlmn(:)
 real(dp) :: cgtmp(2)
 real(dp),parameter :: qphon(3)=(/zero,zero,zero/)
 real(dp),allocatable :: favg_dum(:),fcart_dum(:,:),forold_dum(:,:),fred_tmp(:,:)
 real(dp),allocatable :: gresid_dum(:,:),grhf_dum(:,:),grxc_dum(:,:)
 real(dp),allocatable :: rhog_ep(:,:),scocc(:),str_tmp(:),synlgr_dum(:,:)
 real(dp) :: nhatgr(0,0,0)
 type(pawcprj_type),allocatable :: cprj_tmp(:,:)
 type(pawrhoij_type),allocatable :: pawrhoij_tmp(:)

! *************************************************************************

!Compatibility tests
 if (dtset%positron==0) then
   MSG_BUG('Not valid for dtset%positron=0!')
 end if

 if (istep>1.and.nfft/=electronpositron%nfft) then
   MSG_BUG('Invalid value for nfft!')
 end if

 if (dtset%usewvl==1) then
   MSG_BUG('Not valid for wavelets!')
 end if

 if (dtset%positron==1) then
   do isppol=1,dtset%nsppol
     do ikpt=1,dtset%nkpt
       if (dtset%nband(ikpt+dtset%nkpt*(isppol-1))/=dtset%nband(1)) then
         message = "dtset%positron needs nband to be the same at each k-point !"
         MSG_ERROR(message)
       end if
     end do
   end do
 end if

 comm_cell = mpi_enreg%comm_cell

!-----------------------------------------------------------------------
!Compute new value for calctype (type of electron-positron calculation)
!-----------------------------------------------------------------------
 icalctype0=electronpositron%calctype

 new_calctype=.false.
 if (dtset%positron==1.or.dtset%positron==2) then
   if (istep==1) new_calctype=.true.
 else if (dtset%positron<0) then
   if (ifirst_gs/=0.and.istep==1.and.(.not.always_restart)) new_calctype=.true.
   if (electronpositron%scf_converged) new_calctype=.true.
 end if

!Comment:
!history_level=-1:  not used
!history_level= 0:  rhor from scratch, rhor_ep from scratch or read
!history_level= 1:  rhor in memory, rhor_ep from scratch or read
!history_level= 2:  rhor_ep <-rhor, rhor from scratch
!history_level= 3:  rhor_ep <-> rhor
!history_level= 4:  rhor in memory, rhor_ep in memory
 history_level=-1
 if (dtset%positron==1.or.dtset%positron==2) then
   if (ifirst_gs==0.and.istep==1) history_level=0
   if (ifirst_gs/=0.and.istep==1) history_level=4
 else if (dtset%positron<0) then
   if (.not.electronpositron%scf_converged) then
     if (ifirst_gs/=0.and.istep==1.and.(.not.always_restart)) history_level=4
   else if (electronpositron%scf_converged) then
     if (icalctype0==0) history_level=2
     if (icalctype0> 0) history_level=3
   end if
 end if

 electronpositron%calctype=icalctype0
 if (dtset%positron==1.or.dtset%positron==2) then
   electronpositron%calctype=dtset%positron
 else if (dtset%positron<0) then
   if (electronpositron%scf_converged) then
     if (icalctype0==0) electronpositron%calctype=1
     if (icalctype0>0 ) electronpositron%calctype=3-electronpositron%calctype
   else if (ifirst_gs/=0.and.istep==1) then
     if (always_restart) then
       electronpositron%calctype=0
     else
       electronpositron%calctype=2
!       if (electronpositron%particle==EP_POSITRON) electronpositron%calctype=1
     end if
   end if
 end if

 electronpositron%scf_converged=.false.
 if (new_calctype) electronpositron%istep=electronpositron%istep+1
 if (istep==1) electronpositron%istep=1
 ireadwf=dtfil%ireadwf;if (electronpositron%istep>1) ireadwf=0

!============================================
!The following lines occur only when the type
!of electron-positron calculation changes
!============================================
 if (new_calctype) then

!  Reset some indexes
   if (electronpositron%calctype==0) then
     electronpositron%particle=EP_NOTHING
   else if (electronpositron%calctype==1) then
     electronpositron%particle=EP_ELECTRON
   else if (electronpositron%calctype==2) then
     electronpositron%particle=EP_POSITRON
   end if
   electronpositron%has_pos_ham=mod(electronpositron%calctype,2)
   electronpositron%istep_scf=1
   istep_mix=1

!  -----------------------------------------------------------------------------------------
!  Update forces and stresses
!  If electronpositron%calctype==1: fred_ep/stress_ep are the electronic fred/stress
!  If electronpositron%calctype==2: fred_ep/stress_ep are the positronic fred/stress
!  -----------------------------------------------------------------------------------------
   if (history_level==2.or.history_level==3) then
     optstr=0;optfor=0
     if (allocated(electronpositron%stress_ep)) optstr=stress_needed
     if (allocated(electronpositron%fred_ep).and.forces_needed==2) optfor=1
     if (optfor>0.or.optstr>0) then
       ABI_ALLOCATE(favg_dum,(3))
       ABI_ALLOCATE(fcart_dum,(3,dtset%natom))
       ABI_ALLOCATE(forold_dum,(3,dtset%natom))
       ABI_ALLOCATE(gresid_dum,(3,dtset%natom))
       ABI_ALLOCATE(grhf_dum,(3,dtset%natom))
       ABI_ALLOCATE(grxc_dum,(3,dtset%natom))
       ABI_ALLOCATE(synlgr_dum,(3,dtset%natom))
       ABI_ALLOCATE(fred_tmp,(3,dtset%natom))
       ABI_ALLOCATE(str_tmp,(6))
       forold_dum=zero;n3xccc0=n3xccc
       icalctype=electronpositron%calctype;electronpositron%calctype=-icalctype0
       if (electronpositron%calctype==0) electronpositron%calctype=-100
       if (electronpositron%calctype==-1) n3xccc0=0  ! Note: if calctype=-1, previous calculation was positron
       call forstr(atindx1,cg,cprj,diffor_dum,dtefield,dtset,eigen,electronpositron,energies,&
&       favg_dum,fcart_dum,fock,forold_dum,fred_tmp,grchempottn,gresid_dum,grewtn,grhf_dum,grvdw,grxc_dum,gsqcut,&
&       indsym,kg,kxc,maxfor_dum,mcg,mcprj,mgfft,mpi_enreg,my_natom,n3xccc0,nattyp,nfft,ngfft,&
&       ngrvdw,nhat,nkxc,npwarr,dtset%ntypat,nvresid,occ,optfor,optres,paw_ij,pawang,pawfgr,&
&       pawfgrtab,pawrad,pawrhoij,pawtab,ph1dc,ph1d,psps,rhog,rhor,rprimd,optstr,strsxc,str_tmp,symrec,&
&       synlgr_dum,ucvol,usecprj,vhartr,vpsp,vxc,wvl,xccc3d,xred,ylm,ylmgr,0.0_dp)
       electronpositron%calctype=icalctype
       if (optfor>0) electronpositron%fred_ep(:,:)=fred_tmp(:,:)
       if (optstr>0) electronpositron%stress_ep(:)=str_tmp(:)
       ABI_DEALLOCATE(favg_dum)
       ABI_DEALLOCATE(fcart_dum)
       ABI_DEALLOCATE(forold_dum)
       ABI_DEALLOCATE(gresid_dum)
       ABI_DEALLOCATE(grhf_dum)
       ABI_DEALLOCATE(grxc_dum)
       ABI_DEALLOCATE(synlgr_dum)
       ABI_DEALLOCATE(fred_tmp)
       ABI_DEALLOCATE(str_tmp)
     end if
     if (optfor==0.and.forces_needed>0.and.allocated(electronpositron%fred_ep)) then
       electronpositron%fred_ep(:,:)=fred(:,:)-electronpositron%fred_ep(:,:)
     end if
   end if

!  ----------------------------------------------------------------------------------------------------
!  Initialize/Update densities
!  If electronpositron%calctype==1: rhor is the positronic density, rhor_ep is the electronic density
!  If electronpositron%calctype==2: rhor is the electronic density, rhor_ep is the positronic density
!  ---------------------------------------------------------------------------------------------------
   ABI_ALLOCATE(rhog_ep,(2,nfft))

!  ===== PREVIOUS DENSITY RHOR_EP:
   if (history_level==0.or.history_level==1) then
!    ----- Read from disk
     if (dtset%positron>0) then
       rdwrpaw=dtset%usepaw
       fname=trim(dtfil%fildensin);if (dtset%positron==2) fname=trim(dtfil%fildensin)//'_POSITRON'
       call read_rhor(trim(fname), cplex1, dtset%nspden, nfft, ngfft, rdwrpaw, mpi_enreg, electronpositron%rhor_ep, &
       hdr_den, electronpositron%pawrhoij_ep, comm_cell, check_hdr=hdr)
       etotal_read = hdr_den%etot; call hdr_free(hdr_den)
       call fourdp(1,rhog_ep,electronpositron%rhor_ep,-1,mpi_enreg,nfft,1,ngfft,0)
       if (dtset%usepaw==1.and.allocated(electronpositron%nhat_ep)) then
         call pawmknhat(occtmp,1,0,0,0,0,gprimd,my_natom,dtset%natom,nfft,ngfft,0,&
&         dtset%nspden,dtset%ntypat,pawang,pawfgrtab,nhatgr,electronpositron%nhat_ep,&
&         electronpositron%pawrhoij_ep,electronpositron%pawrhoij_ep,pawtab,&
&         qphon,rprimd,ucvol,dtset%usewvl,xred,distribfft=mpi_enreg%distribfft,&
&         comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&         comm_fft=mpi_enreg%comm_fft,paral_kgb=dtset%paral_kgb,me_g0=mpi_enreg%me_g0)
       end if
     end if
!    ----- Electronic from scratch
     if (dtset%positron<0.and.electronpositron%calctype==1) then
       ecut_eff=dtset%pawecutdg*(dtset%dilatmx)**2
       call getcut(boxcut_dum,ecut_eff,gmet,gsqcut_eff,dtset%iboxcut,std_out,qphon,ngfft)
       call initro(atindx,dtset%densty,gmet,gsqcut_eff,dtset%usepaw,mgfft,mpi_enreg,&
&       psps%mqgrid_vl,dtset%natom,nattyp,nfft,ngfft,dtset%nspden,dtset%ntypat,&
&       psps,pawtab,ph1d,psps%qgrid_vl,rhog_ep,electronpositron%rhor_ep,&
&       dtset%spinat,ucvol,dtset%usepaw,dtset%ziontypat,dtset%znucl)
       if (dtset%usepaw==1) then
         if (size(electronpositron%pawrhoij_ep)>0) then
           ABI_DATATYPE_ALLOCATE(pawrhoij_tmp,(my_natom))
           call initrhoij(electronpositron%pawrhoij_ep(1)%cplex_rhoij,dtset%lexexch,&
&           dtset%lpawu,my_natom,dtset%natom,dtset%nspden,&
&           electronpositron%pawrhoij_ep(1)%nspinor,dtset%nsppol,&
&           dtset%ntypat,pawrhoij_tmp,dtset%pawspnorb,pawtab,cplex1,dtset%spinat,dtset%typat,&
&           ngrhoij=electronpositron%pawrhoij_ep(1)%ngrhoij,&
&           nlmnmix=electronpositron%pawrhoij_ep(1)%lmnmix_sz,&
&           use_rhoij_=electronpositron%pawrhoij_ep(1)%use_rhoij_,&
&           use_rhoijres=electronpositron%pawrhoij_ep(1)%use_rhoijres,&
&           comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
           if (electronpositron%pawrhoij_ep(1)%lmnmix_sz>0) then
             do iatom=1,my_natom
               pawrhoij_tmp(iatom)%kpawmix(:)=electronpositron%pawrhoij_ep(iatom)%kpawmix(:)
             end do
           end if
           call pawrhoij_copy(pawrhoij_tmp,electronpositron%pawrhoij_ep)
           call pawrhoij_free(pawrhoij_tmp)
           ABI_DATATYPE_DEALLOCATE(pawrhoij_tmp)
         end if
         if (allocated(electronpositron%nhat_ep)) then
           call pawmknhat(occtmp,1,0,0,0,0,gprimd,my_natom,dtset%natom,nfft,ngfft,0,&
&           dtset%nspden,dtset%ntypat,pawang,pawfgrtab,nhatgr,electronpositron%nhat_ep,&
&           electronpositron%pawrhoij_ep,electronpositron%pawrhoij_ep,pawtab,qphon,rprimd,&
&           ucvol,dtset%usewvl,xred,distribfft=mpi_enreg%distribfft,&
&           comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&           comm_fft=mpi_enreg%comm_fft,paral_kgb=dtset%paral_kgb,me_g0=mpi_enreg%me_g0)
         end if
       end if
     end if
!    ----- Positronic from scratch
     if (dtset%positron<0.and.electronpositron%calctype==2) then
       electronpositron%rhor_ep(:,1)=one/ucvol
       if (dtset%nspden>=2) electronpositron%rhor_ep(:,2)=half/ucvol
       if (dtset%nspden==4) electronpositron%rhor_ep(:,3:4)=zero
       rhog_ep=zero;rhog_ep(1,1)=one/ucvol
       if (dtset%usepaw==1) then
         do iatom=1,dtset%natom
           electronpositron%pawrhoij_ep(iatom)%rhoijp(:,:)=zero
           electronpositron%pawrhoij_ep(iatom)%nrhoijsel=0
         end do
         if (allocated(electronpositron%nhat_ep)) electronpositron%nhat_ep=zero
       end if
     end if
   end if
!  ----- Deduced from rhor in memory
   if (history_level==2) then
     electronpositron%rhor_ep(:,:)=rhor(:,:)
     rhog_ep(:,:)=rhog(:,:)
     if (dtset%usepaw==1) then
       call pawrhoij_copy(pawrhoij,electronpositron%pawrhoij_ep)
       if (allocated(electronpositron%nhat_ep)) electronpositron%nhat_ep(:,:)=nhat(:,:)
     end if
   end if

!  ===== CURRENT DENSITY RHOR:
   if (history_level==0.or.history_level==2) then
     if (ireadwf==0) then
!      ----- Positronic from scratch
       if (electronpositron%calctype==1) then
         rhor(:,1)=one/ucvol
         if (dtset%nspden>=2) rhor(:,2)=half/ucvol
         if (dtset%nspden==4) rhor(:,3:4)=zero
         rhog=zero;rhog(1,1)=one/ucvol
         if (dtset%usepaw==1) then
           do iatom=1,my_natom
             pawrhoij(iatom)%rhoijp(:,:)=zero
             pawrhoij(iatom)%nrhoijsel=0
           end do
           nhat(:,:)=zero
         end if
       end if
!      ----- Electronic from scratch
       if (electronpositron%calctype==2) then
         ecut_eff=dtset%pawecutdg*(dtset%dilatmx)**2
         call getcut(boxcut_dum,ecut_eff,gmet,gsqcut_eff,dtset%iboxcut,std_out,qphon,ngfft)
         call initro(atindx,dtset%densty,gmet,gsqcut_eff,dtset%usepaw,mgfft,mpi_enreg,&
&         psps%mqgrid_vl,dtset%natom,nattyp,nfft,ngfft,dtset%nspden,dtset%ntypat,&
&         psps,pawtab,ph1d,psps%qgrid_vl,rhog,rhor,dtset%spinat,ucvol,&
&         dtset%usepaw,dtset%ziontypat,dtset%znucl)

         if (dtset%usepaw==1) then
           if (size(pawrhoij)>0) then
             ABI_DATATYPE_ALLOCATE(pawrhoij_tmp,(my_natom))
             call initrhoij(pawrhoij(1)%cplex_rhoij,dtset%lexexch,dtset%lpawu,&
&             my_natom,dtset%natom,dtset%nspden,pawrhoij(1)%nspinor,dtset%nsppol,&
&             dtset%ntypat,pawrhoij_tmp,dtset%pawspnorb,pawtab,pawrhoij(1)%qphase,dtset%spinat,&
&             dtset%typat,ngrhoij=pawrhoij(1)%ngrhoij,nlmnmix=pawrhoij(1)%lmnmix_sz,&
&             use_rhoij_=pawrhoij(1)%use_rhoij_,use_rhoijres=pawrhoij(1)%use_rhoijres,&
&             comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
             do iatom=1,my_natom
               pawrhoij_tmp(iatom)%kpawmix(:)=pawrhoij(iatom)%kpawmix(:)
             end do
             call pawrhoij_copy(pawrhoij_tmp,pawrhoij)
             call pawrhoij_free(pawrhoij_tmp)
             ABI_DATATYPE_DEALLOCATE(pawrhoij_tmp)
           end if
           call pawmknhat(occtmp,1,0,0,0,0,gprimd,my_natom,dtset%natom,nfft,ngfft,0,&
&           dtset%nspden,dtset%ntypat,pawang,pawfgrtab,nhatgr,nhat,&
&           pawrhoij,pawrhoij,pawtab,qphon,rprimd,ucvol,dtset%usewvl,xred, &
&           comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&           comm_fft=mpi_enreg%comm_fft,paral_kgb=dtset%paral_kgb,&
&           me_g0=mpi_enreg%me_g0,distribfft=mpi_enreg%distribfft)
         end if

       end if
     end if
   end if

!  ===== EXCHANGE POSITRONIC AND ELECTRONIC DENSITY (CURRENT AND PREVIOUS)
   if (history_level==3) then
     do ispden=1,dtset%nspden
       do ifft=1,nfft
         rhotmp=rhor(ifft,ispden)
         rhor(ifft,ispden)=electronpositron%rhor_ep(ifft,ispden)
         electronpositron%rhor_ep(ifft,ispden)=rhotmp
       end do
     end do
     rhog_ep(:,:)=rhog
     call fourdp(1,rhog,rhor,-1,mpi_enreg,nfft,1,ngfft,0)
!    If PAW, exchange "positronic" and "electronic" rhoij
     if (dtset%usepaw==1) then
       if (size(pawrhoij)>0.and.size(electronpositron%pawrhoij_ep)>0) then
         ABI_DATATYPE_ALLOCATE(pawrhoij_tmp,(my_natom))
         call pawrhoij_alloc(pawrhoij_tmp,pawrhoij(1)%cplex_rhoij,pawrhoij(1)%nspden,&
&         pawrhoij(1)%nspinor,pawrhoij(1)%nsppol,dtset%typat,&
&         pawtab=pawtab,ngrhoij=pawrhoij(1)%ngrhoij,nlmnmix=pawrhoij(1)%lmnmix_sz,&
&         use_rhoij_=pawrhoij(1)%use_rhoij_,use_rhoijres=pawrhoij(1)%use_rhoijres, &
&         comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
         call pawrhoij_copy(pawrhoij,pawrhoij_tmp)
         call pawrhoij_copy(electronpositron%pawrhoij_ep,pawrhoij)
         call pawrhoij_copy(pawrhoij_tmp,electronpositron%pawrhoij_ep)
         call pawrhoij_free(pawrhoij_tmp)
         ABI_DATATYPE_DEALLOCATE(pawrhoij_tmp)
       end if
       if (allocated(electronpositron%nhat_ep)) then
         do ispden=1,dtset%nspden
           do ifft=1,nfft
             rhotmp=nhat(ifft,ispden)
             nhat(ifft,ispden)=electronpositron%nhat_ep(ifft,ispden)
             electronpositron%nhat_ep(ifft,ispden)=rhotmp
           end do
         end do
       end if
     end if
   end if

!  ===== COMPUTE HARTREE POTENTIAL ASSOCIATED TO RHOR_EP
   if (history_level==4) then
     call fourdp(1,rhog_ep,electronpositron%rhor_ep,-1,mpi_enreg,nfft,1,ngfft,0)
   end if
   if (history_level/=-1) then
     call hartre(1,gsqcut,dtset%usepaw,mpi_enreg,nfft,ngfft,rhog_ep,rprimd,&
&     electronpositron%vha_ep)
     electronpositron%vha_ep=-electronpositron%vha_ep
   else
     electronpositron%vha_ep=zero
   end if
   ABI_DEALLOCATE(rhog_ep)

!  ----------------------------------------------------------------------
!  Initialize/Update energies
!  ----------------------------------------------------------------------
   electronpositron%etotal_prev=etotal
   electronpositron%maxfor_prev=maxfor

!  Inits/exchange news energies
!  Retrieve energy of non-evolving particle(s)
   if (history_level== 0) then
     call energies_init(energies)
     call energies_init(electronpositron%energies_ep)
     if (dtset%positron>0) energies%e0_electronpositron=etotal_read
     if (dtset%positron<0) energies%e0_electronpositron=zero
   else if (history_level== 1) then
     call energies_init(electronpositron%energies_ep)
     if (dtset%positron>0) energies%e0_electronpositron=etotal_read
   else if (history_level== 2) then
     call energies_copy(energies,electronpositron%energies_ep)
     call energies_init(energies)
     energies%e0_electronpositron=electronpositron%e0
   else if (history_level== 3) then
     call energies_copy(electronpositron%energies_ep,energies_tmp)
     call energies_copy(energies,electronpositron%energies_ep)
     call energies_copy(energies_tmp,energies)
     energies%e0_electronpositron=electronpositron%e0
!    else if (history_level== 4) then
   end if

!  Adjust core psps energy
   if (electronpositron%calctype/=1) energies%e_corepsp=ecore/ucvol

!  -----------------------------------------------------------------------------------------
!  Update wavefunctions
!  If electronpositron%calctype==1: cg are the positronic WFs, cg_ep are the electronic WFs
!  If electronpositron%calctype==2: cg are the electronic WFs, cg_ep are the positronic WFs
!  -----------------------------------------------------------------------------------------
   if (electronpositron%dimcg>0.or.electronpositron%dimcprj>0) then

     if (history_level==0.or.history_level==1) then
       electronpositron%cg_ep=zero
     end if

     if (history_level==2) then
       if (electronpositron%dimcg>0) then
         do icg=1,electronpositron%dimcg
           electronpositron%cg_ep(1:2,icg)=cg(1:2,icg)
         end do
       end if
       if (dtset%usepaw==1.and.electronpositron%dimcprj>0) then
         call pawcprj_copy(cprj,electronpositron%cprj_ep)
       end if
     end if

     if (history_level==3) then
       if (electronpositron%dimcg>0) then
         do icg=1,electronpositron%dimcg
           cgtmp(1:2)=electronpositron%cg_ep(1:2,icg)
           electronpositron%cg_ep(1:2,icg)=cg(1:2,icg)
           cg(1:2,icg)=cgtmp(1:2)
         end do
       end if
       if (dtset%usepaw==1.and.electronpositron%dimcprj>0) then
         ABI_ALLOCATE(nlmn,(dtset%natom))
         ABI_DATATYPE_ALLOCATE(cprj_tmp,(dtset%natom,electronpositron%dimcprj))
         do iatom=1,dtset%natom
           nlmn(iatom)=cprj(iatom,1)%nlmn
         end do
         call pawcprj_alloc(cprj_tmp,cprj(1,1)%ncpgr,nlmn)
         ABI_DEALLOCATE(nlmn)
         call pawcprj_copy(electronpositron%cprj_ep,cprj_tmp)
         call pawcprj_copy(cprj,electronpositron%cprj_ep)
         call pawcprj_copy(cprj_tmp,cprj)
         call pawcprj_free(cprj_tmp)
         ABI_DATATYPE_DEALLOCATE(cprj_tmp)
       end if
     end if

   end if ! dimcg>0 or dimcprj>0

!  -----------------------------------------------------------------------------------------------------------
!  Initialize/Update occupations
!  If electronpositron%calctype==1: occ are the positronic occupations, occ_ep are the electronic occupations
!  If electronpositron%calctype==2: occ are the electronic occupations, occ_ep are the positronic occupations
!  -----------------------------------------------------------------------------------------------------------
!  When needed, precompute electronic occupations with semiconductor occupancies
   need_scocc=.false.
   if (electronpositron%dimocc>0.and.electronpositron%calctype==1.and. &
&   (history_level==0.or.history_level==1)) need_scocc=.true.
   if (electronpositron%calctype==2.and.ireadwf==0.and. &
&   (history_level==0.or.history_level==2.or. &
&   (history_level==3.and.electronpositron%dimocc==0))) need_scocc=.true.
   if (need_scocc) then
     nelect=-dtset%charge
     do iatom=1,dtset%natom
       nelect=nelect+dtset%ziontypat(dtset%typat(iatom))
     end do
     maxocc=two/real(dtset%nsppol*dtset%nspinor,dp)
     nocc=(nelect-tol8)/maxocc + 1
     nocc=min(nocc,dtset%nband(1)*dtset%nsppol)
     occlast=nelect-maxocc*(nocc-1)
     ABI_ALLOCATE(scocc,(dtset%nband(1)*dtset%nsppol))
     scocc=zero
     if (1<nocc)  scocc(1:nocc-1)=maxocc
     if (1<=nocc) scocc(nocc)=occlast
   end if

!  ===== PREVIOUS OCCUPATIONS OCC_EP:
   if (electronpositron%dimocc>0) then
     if (history_level==0.or.history_level==1) then
!      ----- Electronic from scratch
       if (electronpositron%calctype==1) then
!        Initialize electronic occupations with semiconductor occupancies
         do ikpt=1,dtset%nkpt
           do iband=1,dtset%nband(1)
             do isppol=1,dtset%nsppol
               electronpositron%occ_ep(iband+dtset%nband(1)*(ikpt-1+dtset%nkpt*(isppol-1)))=&
&               scocc(isppol+dtset%nsppol*(iband-1))
             end do
           end do
         end do
       end if
!      ----- Positronic from scratch
       if (electronpositron%calctype==1) then
!        Initialize positronic occupations with only one positron (or less)
         electronpositron%occ_ep(:)=zero
         isppol=1;iocc=1
         do ikpt=1,dtset%nkpt
           electronpositron%occ_ep(iocc)=electronpositron%posocc
           iocc=iocc+dtset%nband(ikpt+dtset%nkpt*(isppol-1))
         end do
       end if
     end if
!    ----- Deduced from occ in memory
     if (history_level==2) then
       electronpositron%occ_ep(:)=occ(:)
     end if
   end if ! dimocc>0

!  ===== CURRENT OCCUPATIONS OCC:
   if (history_level==0.or.history_level==2.or.(history_level==3.and.electronpositron%dimocc==0)) then
     if (ireadwf==0) then
!      ----- Positronic from scratch
       if (electronpositron%calctype==1) then
!        Initialize positronic occupations with only one positron (or less)
         occ(:)=zero
         isppol=1;iocc=1
         do ikpt=1,dtset%nkpt
           occ(iocc)=electronpositron%posocc
           iocc=iocc+dtset%nband(ikpt+dtset%nkpt*(isppol-1))
         end do
       end if
!      ----- Electronic from scratch
       if (electronpositron%calctype==2) then
!        Initialize electronic occupations with semiconductor occupancies
         do ikpt=1,dtset%nkpt
           do iband=1,dtset%nband(1)
             do isppol=1,dtset%nsppol
               occ(iband+dtset%nband(1)*(ikpt-1+dtset%nkpt*(isppol-1)))=&
&               scocc(isppol+dtset%nsppol*(iband-1))
             end do
           end do
         end do
       end if
     end if
   end if

!  ===== EXCHANGE POSITRONIC AND ELECTRONIC OCCUPATIONS (CURRENT AND PREVIOUS)
   if (history_level==3.and.electronpositron%dimocc>0) then
     do iocc=1,electronpositron%dimocc
       occtmp=occ(iocc)
       occ(iocc)=electronpositron%occ_ep(iocc)
       electronpositron%occ_ep(iocc)=occtmp
     end do
   end if

   if (need_scocc)  then
     ABI_DEALLOCATE(scocc)
   end if

!  -----------------------------------------------------------------------------------------------------------
!  Initialize/Update eigen energies
!  If electronpositron%calctype==1: eigen are the positronic eigen E, eigen_ep are the electronic eigen E
!  If electronpositron%calctype==2: eigen are the electronic eigen E, eigen_ep are the positronic eigen E
!  -----------------------------------------------------------------------------------------------------------

!  ===== PREVIOUS EIGEN ENERGIES EIGEN_EP:
   if (electronpositron%dimeigen>0) then
     if (history_level==0.or.history_level==1) then
!      ----- Electronic or positronic from scratch
       electronpositron%eigen_ep(:)=zero
     end if
!    ----- Deduced from eigen in memory
     if (history_level==2) then
       electronpositron%eigen_ep(:)=eigen(:)
     end if
   end if ! dimeigen>0

!  ===== CURRENT EIGEN ENERGIES EIGEN:
   if (history_level==0.or.history_level==2.or.(history_level==3.and.electronpositron%dimeigen==0)) then
     if (ireadwf==0) then
!      ----- Electronic or positronic from scratch
       eigen(:)=zero
     end if
   end if

!  ===== EXCHANGE POSITRONIC AND ELECTRONIC EIGEN ENERGIES (CURRENT AND PREVIOUS)
   if (history_level==3.and.electronpositron%dimeigen>0) then
     do iocc=1,electronpositron%dimeigen
       eigtmp=eigen(iocc)
       eigen(iocc)=electronpositron%eigen_ep(iocc)
       electronpositron%eigen_ep(iocc)=eigtmp
     end do
   end if

!  =============================================
 end if  ! the type of e-p calculation changes
!=============================================

!------------------------------------------------------------------
!Write messages
!------------------------------------------------------------------
 if (istep_mix==1.and.dtset%positron/=0) then
!  Log message
   if (electronpositron%calctype==0) then
     message = 'Were are now performing an electronic ground-state calculation...'
   else if (electronpositron%calctype==1) then
     message = 'Were are now performing a positronic ground-state calculation...'
   else if (electronpositron%calctype==2) then
     message = 'Were are now performing an electronic ground-state calculation in presence of a positron...'
   end if
   MSG_COMMENT(message)
!  Output message
   if (dtset%positron<0) then
     if (electronpositron%calctype==0) then
       TypeCalcStrg='ELECTRONIC GROUND-STATE CALCULATION'
     else if (electronpositron%calctype==1) then
       TypeCalcStrg='POSITRONIC GROUND-STATE CALCULATION IN PRESENCE OF ELECTRONS AND IONS'
     else if (electronpositron%calctype==2) then
       TypeCalcStrg='ELECTRONIC GROUND-STATE CALCULATION IN PRESENCE OF A POSITRON'
     end if
     if (istep>1) then
       write(message,'(2a,i3,2a)') ch10,'TC-DFT STEP ',electronpositron%istep,' - ',trim(TypeCalcStrg)
     else
       write(message,'(a,i3,2a)') 'TC-DFT STEP ',electronpositron%istep,' - ',trim(TypeCalcStrg)
     end if
     call wrtout(ab_out,message,'COLL')
   end if
 end if

end subroutine setup_positron
!!***

!!****f* ABINIT/poslifetime
!! NAME
!! poslifetime
!!
!! FUNCTION
!! Calculate the positron lifetime
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | nspden=number of spin-density components
!!   | ntypat=number of atom types
!!   | paral_kgb=flag controlling (k,g,bands) parallelization
!!   | pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!   | usepaw=flag for PAW
!!  gprimd(3,3)= dimensional reciprocal space primitive translations
!!  mpi_enreg= information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  n3xccc= dimension of the xccc3d array (0 or nfft).
!!  nfft= number of FFT grid points
!!  ngfft(18)= contain all needed information about 3D FFT
!!  nhat(nfft,nspden)=charge compensation density (content depends on electronpositron%particle)
!!  option= if 1, calculate positron lifetime for whole density
!!          if 2, calculate positron lifetime for given state
!!          if 3, calculate positron lifetime for given state with IPM
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(my_natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  rhor(nfft,nspden)=total electron/positron density (content depends on electronpositron%particle)
!!  ucvol=unit cell volume in bohr**3.
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  ===== Optional arguments, used only if option>1 =====
!!  pawrhoij_dop_el(my_natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies of one state
!!  rhor_dop_el(nfft)=electron density of given state for the state dependent scheme
!!  ===== Optional argument =====
!!  pawrhoij_ep(my_natom*usepaw) <type(pawrhoij_type)>= atomic occupancies to be used in place of
!!                                                      electronpositron%pawrhoij_ep
!!
!! OUTPUT
!!  rate= annihilation rate of a given state needed for state dependent scheme for doppler broadening
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!
!! PARENTS
!!      outscfcv,posdoppler
!!
!! CHILDREN
!!      gammapositron,gammapositron_fft,mkdenpos,nderiv_gen,pawdensities
!!      pawxcsum,simp_gen,wrtout,xmpi_sum
!!
!! SOURCE

subroutine poslifetime(dtset,electronpositron,gprimd,my_natom,mpi_enreg,n3xccc,nfft,ngfft,nhat,&
&                      option,pawang,pawrad,pawrhoij,pawtab,rate,rate_paw,rhor,ucvol,xccc3d,&
&                      rhor_dop_el,pawrhoij_dop_el,pawrhoij_ep) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_natom,n3xccc,nfft,option
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: rate,rate_paw
 type(dataset_type), intent(in) :: dtset
 type(electronpositron_type),pointer :: electronpositron
 type(MPI_type),intent(in) :: mpi_enreg
 type(pawang_type), intent(in) :: pawang
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3),nhat(nfft,dtset%nspden*dtset%usepaw),xccc3d(n3xccc)
 real(dp),intent(in),target :: rhor(nfft,dtset%nspden)
 real(dp),optional,intent(in) :: rhor_dop_el(nfft)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*dtset%usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom*dtset%usepaw)
 type(pawrhoij_type),optional,intent(in) :: pawrhoij_dop_el(my_natom*dtset%usepaw)
 type(pawrhoij_type),optional,target,intent(in) :: pawrhoij_ep(my_natom*dtset%usepaw)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer :: cplex,iatom,ierr,ifft,igam,ii,ilm,ilm1,ilm2,iloop,ipt,ir,isel
 integer :: itypat,iwarn,iwarnj,iwarnp,lm_size,lmn2_size,mesh_size
 integer :: nfftot,ngamma,ngr,ngrad,nspden_ep,opt_dens,usecore
 logical,parameter :: include_nhat_in_gamma=.false.
 real(dp),parameter :: delta=1.d-4
 real(dp) :: fact,fact2,intg
 real(dp) :: lambda_core    ,lambda_core_ipm    ,lambda    ,lambda_ipm
 real(dp) :: lambda_core_paw,lambda_core_paw_ipm,lambda_paw,lambda_paw_ipm
 real(dp) :: lifetime,lifetime_ipm,nbec,nbev,nbp,rdum,sqfpi,units
 character(len=500) :: msg
!arrays
 integer,allocatable :: igamma(:)
 logical,allocatable :: lmselect(:),lmselect_ep(:),lmselect_dum(:)
 real(dp) :: mpibuf(4)
 real(dp),parameter :: qphon(3)=(/zero,zero,zero/),lsign(2)=(/one,-one/)
 real(dp),allocatable :: d1gam(:,:,:),d2gam(:,:,:),ff(:),gam_(:,:,:),gamma(:,:),gammam(:,:,:),gg(:,:)
 real(dp),allocatable :: grhocore2(:),grhocor2_(:),grhoe2(:),grho2_(:)
 real(dp),allocatable :: nhat1(:,:,:),nhat1_ep(:,:,:),nhat1_j(:,:,:)
 real(dp),allocatable :: rho_(:),rho_ep_(:),rho1(:,:,:),rho1_ep(:,:,:),rho1_j(:,:,:)
 real(dp),allocatable :: rhoarr1(:),rhoarr1_ep(:),rhoarr1_j(:),rhoarr2(:)
 real(dp),allocatable :: rhocore(:),rhocor_(:),rhoe(:,:),rhop(:,:),rhor_dop_el_(:)
 real(dp),allocatable :: rhosph(:),rhosph_ep(:),rhosph_j(:),rhotot(:,:),rhotot_ep(:,:)
 real(dp),allocatable :: rhotot_j(:,:),trho1(:,:,:),trho1_ep(:,:,:),trho1_j(:,:,:)
 real(dp),allocatable :: v1sum(:,:),v2sum(:,:,:)
 real(dp),pointer :: rhor_(:,:),rhor_ep_(:,:)
 type(pawrhoij_type),pointer :: pawrhoij_ep_(:)

! *************************************************************************

 DBG_ENTER("COLL")

!Tests for developers
 if (.not.associated(electronpositron)) then
   msg='electronpositron variable must be associated!'
   MSG_BUG(msg)
 end if
 if (option/=1) then
   if ((.not.present(rhor_dop_el)).or.(.not.present(pawrhoij_dop_el))) then
     msg='when option/=1, rhor_dop_el and pawrhoij_dop_el must be present!'
     MSG_BUG(msg)
   end if
 end if

!Constants
 fact=0.0
 cplex=1;nspden_ep=1
 usecore=n3xccc/nfft
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 ngrad=1;if (electronpositron%ixcpositron==3.or.electronpositron%ixcpositron==31) ngrad=2
 iwarn=0;iwarnj=0;iwarnp=1
 sqfpi=sqrt(four_pi)

!Compatibility tests
 if (electronpositron%particle==EP_NOTHING) then
   msg='Not valid for electronpositron%particle=NOTHING!'
   MSG_BUG(msg)
 end if
 if (electronpositron%nfft/=nfft) then
   msg='nfft/=electronpositron%nfft!'
   MSG_BUG(msg)
 end if
 if (dtset%usepaw==1) then
   if(dtset%pawxcdev==0.and.ngrad==2) then
     msg='GGA is not implemented for pawxcdev=0 (use dtset%pawxcdev/=0)!'
     MSG_BUG(msg)
   end if
 end if

!Select type(s) of enhancement factor
 if ((electronpositron%ixcpositron==1.or.electronpositron%ixcpositron==3).and.option==1) then
   ngamma=2
   ABI_ALLOCATE(igamma,(ngamma))
   igamma(1)=1;igamma(2)=2
 else
   ngamma=1
   ABI_ALLOCATE(igamma,(ngamma))
   if (electronpositron%ixcpositron==-1) igamma(1)=0
   if (electronpositron%ixcpositron== 2) igamma(1)=4
   if (electronpositron%ixcpositron==11.or.electronpositron%ixcpositron==31) igamma(1)=3
   if (electronpositron%ixcpositron==1.or.electronpositron%ixcpositron==3) igamma(1)=2
 end if

!Select density according to nhat choice
 if (dtset%usepaw==0.or.include_nhat_in_gamma) then
   rhor_ => rhor
   rhor_ep_ => electronpositron%rhor_ep
 else
   ABI_ALLOCATE(rhor_,(nfft,dtset%nspden))
   ABI_ALLOCATE(rhor_ep_,(nfft,dtset%nspden))
   rhor_=rhor-nhat
   rhor_ep_=electronpositron%rhor_ep-electronpositron%nhat_ep
 end if

!Eventually overwrite electronpositron%pawrhoij_ep
 if (present(pawrhoij_ep)) then
   pawrhoij_ep_ => pawrhoij_ep
 else
   pawrhoij_ep_ => electronpositron%pawrhoij_ep
 end if

!Loop on different enhancement factors
 do igam=1,ngamma

!  Compute electron-positron annihilation rate using pseudo densities (plane waves)
!  ----------------------------------------------------------------------------------------

!  Select the densities and make them positive
   ABI_ALLOCATE(rhoe,(nfft,nspden_ep))
   ABI_ALLOCATE(rhop,(nfft,nspden_ep))
   if (electronpositron%particle==EP_ELECTRON) then
     rhoe(:,1)=rhor_ep_(:,1);rhop(:,1)=rhor_(:,1)
   else if (electronpositron%particle==EP_POSITRON) then
     rhoe(:,1)=rhor_(:,1);rhop(:,1)=rhor_ep_(:,1)
   end if
   call mkdenpos(iwarn ,nfft,nspden_ep,1,rhoe,dtset%xc_denpos)
   call mkdenpos(iwarnp,nfft,nspden_ep,1,rhop,dtset%xc_denpos)
   if (option/=1) then
     ABI_ALLOCATE(rhor_dop_el_,(nfft))
     rhor_dop_el_(:)=rhor_dop_el(:)
     call mkdenpos(iwarnp,nfft,1,1,rhor_dop_el_,dtset%xc_denpos)
   end if

!  Compute enhancement factor at each FFT grid point
!  gamma(:,1): using total   electronic density
!  gamma(:,2): using valence electronic density
   ABI_ALLOCATE(gamma,(nfft,2))
   if (option==1.or.option==2) then
     call gammapositron_fft(electronpositron,gamma,gprimd,igamma(igam),mpi_enreg,&
&     n3xccc,nfft,ngfft,rhoe,rhop,xccc3d)
   else
     gamma=one
   end if

!  Compute positron annihilation rates
   lambda     =zero;lambda_ipm     =zero
   lambda_core=zero;lambda_core_ipm=zero
   if (option==1) then
     do ifft=1,nfft
       lambda    =lambda    +rhop(ifft,1)*rhoe(ifft,1)*gamma(ifft,1)
       lambda_ipm=lambda_ipm+rhop(ifft,1)*rhoe(ifft,1)*gamma(ifft,2)
     end do
   else
     do ifft=1,nfft
       lambda    =lambda    +rhop(ifft,1)*rhor_dop_el_(ifft)*gamma(ifft,1)
       lambda_ipm=lambda_ipm+rhop(ifft,1)*rhor_dop_el_(ifft)*gamma(ifft,2)
     end do
   end if
   if (usecore==1) then
     do ifft=1,nfft
       lambda_core    =lambda_core    +rhop(ifft,1)*xccc3d(ifft)*gamma(ifft,1)
       lambda_core_ipm=lambda_core_ipm+rhop(ifft,1)*xccc3d(ifft)
     end do
   end if
   lambda         =lambda         *ucvol/dble(nfftot)
   lambda_ipm     =lambda_ipm     *ucvol/dble(nfftot)
   lambda_core    =lambda_core    *ucvol/dble(nfftot)
   lambda_core_ipm=lambda_core_ipm*ucvol/dble(nfftot)
   ABI_DEALLOCATE(gamma)
   ABI_DEALLOCATE(rhoe)
   ABI_DEALLOCATE(rhop)
   if (option/=1) then
     ABI_DEALLOCATE(rhor_dop_el_)
   end if
!  NC pseudopotential: check electrons/positron number
   if (dtset%usepaw==0.and.igam==ngamma) then
     nbec=zero;nbev=zero;nbp=zero
     if (electronpositron%particle==EP_ELECTRON) then
       do ifft=1,nfft
         nbec=nbec+xccc3d(ifft)
         nbev=nbev+electronpositron%rhor_ep(ifft,1)
         nbp =nbp +rhor(ifft,1)
       end do
     else
       do ifft=1,nfft
         nbec=nbec+xccc3d(ifft)
         nbev=nbev+rhor(ifft,1)
         nbp =nbp +electronpositron%rhor_ep(ifft,1)
       end do
     end if
     nbec=nbec*ucvol/dble(nfftot)
     nbev=nbev*ucvol/dble(nfftot)
     nbp =nbp *ucvol/dble(nfftot)
   end if

!  MPI parallelization
   if(mpi_enreg%nproc_fft>1)then
     call xmpi_sum(lambda    ,mpi_enreg%comm_fft,ierr)
     call xmpi_sum(lambda_ipm,mpi_enreg%comm_fft,ierr)
     call xmpi_sum(lambda_core    ,mpi_enreg%comm_fft,ierr)
     call xmpi_sum(lambda_core_ipm,mpi_enreg%comm_fft,ierr)
     if (dtset%usepaw==0.and.igam==ngamma) then
       call xmpi_sum(nbec,mpi_enreg%comm_fft,ierr)
       call xmpi_sum(nbev,mpi_enreg%comm_fft,ierr)
       call xmpi_sum(nbp ,mpi_enreg%comm_fft,ierr)
     end if
   end if


!  PAW: add on-site contributions to electron-positron annihilation rate
!  ----------------------------------------------------------------------------------------
   if (dtset%usepaw==1) then

     lambda_paw     =zero;lambda_paw_ipm     =zero
     lambda_core_paw=zero;lambda_core_paw_ipm=zero

!    Loop on atoms
     do iatom=1,my_natom

       itypat=pawrhoij(iatom)%itypat
       lmn2_size=pawtab(itypat)%lmn2_size
       mesh_size=pawtab(itypat)%mesh_size
       lm_size=pawtab(itypat)%lcut_size**2
       cplex=1
       ngr=0;if (ngrad==2) ngr=mesh_size

!      Allocations of "on-site" densities
       ABI_ALLOCATE(rho1 ,(cplex*mesh_size,lm_size,nspden_ep))
       ABI_ALLOCATE(trho1,(cplex*mesh_size,lm_size,nspden_ep))
       ABI_ALLOCATE(rho1_ep ,(cplex*mesh_size,lm_size,nspden_ep))
       ABI_ALLOCATE(trho1_ep,(cplex*mesh_size,lm_size,nspden_ep))
       if (option/=1) then
         ABI_ALLOCATE(rho1_j ,(cplex*mesh_size,lm_size,nspden_ep))
         ABI_ALLOCATE(trho1_j,(cplex*mesh_size,lm_size,nspden_ep))
       else
         ABI_ALLOCATE(rho1_j ,(0,0,0))
         ABI_ALLOCATE(trho1_j,(0,0,0))
       end if
       if (include_nhat_in_gamma) then
         ABI_ALLOCATE(nhat1,(cplex*mesh_size,lm_size,nspden_ep))
         ABI_ALLOCATE(nhat1_ep,(cplex*mesh_size,lm_size,nspden_ep))
       else
         ABI_ALLOCATE(nhat1,(0,0,0))
         ABI_ALLOCATE(nhat1_ep,(0,0,0))
       end if
       if (include_nhat_in_gamma.and.option/=1) then
         ABI_ALLOCATE(nhat1_j,(cplex*mesh_size,lm_size,nspden_ep))
       else
         ABI_ALLOCATE(nhat1_j,(0,0,0))
       end if
       ABI_ALLOCATE(lmselect,(lm_size))
       ABI_ALLOCATE(lmselect_ep,(lm_size))
       ABI_ALLOCATE(lmselect_dum,(lm_size))

!      Compute "on-site" densities (n1, ntild1, nhat1) for electron and positron =====
       lmselect(:)=.true.
       opt_dens=1;if (include_nhat_in_gamma) opt_dens=0
       call pawdensities(rdum,cplex,iatom,lmselect,lmselect_dum,lm_size,nhat1,nspden_ep,1,&
&       0,opt_dens,-1,0,pawang,0,pawrad(itypat),pawrhoij(iatom),&
&       pawtab(itypat),rho1,trho1)
       lmselect_ep(:)=.true.
       call pawdensities(rdum,cplex,iatom,lmselect_ep,lmselect_dum,lm_size,nhat1_ep,nspden_ep,1,&
&       0,opt_dens,-1,0,pawang,0,pawrad(itypat),pawrhoij_ep_(iatom),&
&       pawtab(itypat),rho1_ep,trho1_ep)

!      For state dependent scheme in Doppler                                       =====
!      Compute "on-site" densities (n1, ntild1, nhat1) for a given electron state j=====
       if (option/=1) then
         opt_dens=1;if (include_nhat_in_gamma) opt_dens=0
         call pawdensities(rdum,cplex,iatom,lmselect,lmselect_dum,lm_size,nhat1_j,nspden_ep,1,&
&         0,opt_dens,-1,0,pawang,0,pawrad(itypat),pawrhoij_dop_el(iatom),&
&         pawtab(itypat),rho1_j,trho1_j)
       end if

!      Compute contribution to annihilation rate:
!      Loop: first step: compute all-electron contribution (from n^1, n_c)
!      2nd   step: compute pseudo contribution (from tild_n^1, hat_n^1, tild_n_c)
       do iloop=1,2
         if (iloop==1) usecore=1
         if (iloop==2) usecore=pawtab(itypat)%usetcore
         ABI_ALLOCATE(rhocore,(mesh_size))

!        First formalism: use densities on r,theta,phi
         if (dtset%pawxcdev==0) then

           ABI_ALLOCATE(gamma,(mesh_size,2))
           ABI_ALLOCATE(rhoarr1,(mesh_size))
           ABI_ALLOCATE(rhoarr1_ep,(mesh_size))
           if (option/=1) then
             ABI_ALLOCATE(rhoarr1_j,(mesh_size))
           end if
!          Loop on the angular part
           do ipt=1,pawang%angl_size
!            Build densities
             rhoarr1=zero;rhoarr1_ep=zero;rhocore=zero
             if (option/=1) rhoarr1_j=zero
             if (iloop==1) then
               do ilm=1,lm_size
                 if (lmselect(ilm)) rhoarr1(:)=rhoarr1(:)+rho1(:,ilm,1)*pawang%ylmr(ilm,ipt)
               end do
               if (option/=1) then
                 do ilm=1,lm_size
                   if (lmselect(ilm)) rhoarr1_j(:)=rhoarr1_j(:)+rho1_j(:,ilm,1)*pawang%ylmr(ilm,ipt)
                 end do
               end if
               do ilm=1,lm_size
                 if (lmselect_ep(ilm)) rhoarr1_ep(:)=rhoarr1_ep(:)+rho1_ep(:,ilm,1)*pawang%ylmr(ilm,ipt)
               end do
               if (usecore==1) rhocore(:)=pawtab(itypat)%coredens(:)
             else
               if (include_nhat_in_gamma) then
                 do ilm=1,lm_size
                   if (lmselect(ilm)) rhoarr1(:)=rhoarr1(:)+(trho1(:,ilm,1)+nhat1(:,ilm,1))*pawang%ylmr(ilm,ipt)
                 end do
                 if (option/=1) then
                   do ilm=1,lm_size
                     if (lmselect(ilm)) rhoarr1_j(:)=rhoarr1_j(:)+(trho1_j(:,ilm,1)+nhat1_j(:,ilm,1))*pawang%ylmr(ilm,ipt)
                   end do
                 end if
                 do ilm=1,lm_size
                   if (lmselect_ep(ilm)) rhoarr1_ep(:)=rhoarr1_ep(:)+(trho1_ep(:,ilm,1)+nhat1_ep(:,ilm,1))*pawang%ylmr(ilm,ipt)
                 end do
               else
                 do ilm=1,lm_size
                   if (lmselect(ilm)) rhoarr1(:)=rhoarr1(:)+trho1(:,ilm,1)*pawang%ylmr(ilm,ipt)
                 end do
                 if (option/=1) then
                   do ilm=1,lm_size
                     if (lmselect(ilm)) rhoarr1_j(:)=rhoarr1(:)+trho1_j(:,ilm,1)*pawang%ylmr(ilm,ipt)
                   end do
                 end if
                 do ilm=1,lm_size
                   if (lmselect_ep(ilm)) rhoarr1_ep(:)=rhoarr1_ep(:)+trho1_ep(:,ilm,1)*pawang%ylmr(ilm,ipt)
                 end do
               end if
               if (usecore==1) rhocore(:)=pawtab(itypat)%tcoredens(:,1)
             end if
!            Make the densities positive
             if (electronpositron%particle==EP_ELECTRON) then
               call mkdenpos(iwarnp,mesh_size,1,1,rhoarr1   ,dtset%xc_denpos)
               call mkdenpos(iwarn ,mesh_size,1,1,rhoarr1_ep,dtset%xc_denpos)
             else if (electronpositron%particle==EP_POSITRON) then
               call mkdenpos(iwarn ,mesh_size,1,1,rhoarr1   ,dtset%xc_denpos)
               call mkdenpos(iwarnp,mesh_size,1,1,rhoarr1_ep,dtset%xc_denpos)
               if (option/=1) then
                 call mkdenpos(iwarnj,mesh_size,1,1,rhoarr1_j,dtset%xc_denpos)
               end if
             end if
!            Compute Gamma
             ABI_ALLOCATE(grhoe2,(ngr))
             ABI_ALLOCATE(grhocore2,(ngr))
             if (option==1.or.option==2) then
               if (electronpositron%particle==EP_ELECTRON) then
                 call gammapositron(gamma,grhocore2,grhoe2,igamma(igam),ngr,mesh_size,&
&                 rhocore,rhoarr1_ep,rhoarr1,usecore)
               else if (electronpositron%particle==EP_POSITRON) then
                 call gammapositron(gamma,grhocore2,grhoe2,igamma(igam),ngr,mesh_size,&
&                 rhocore,rhoarr1,rhoarr1_ep,usecore)
               end if
             else
               gamma(:,:)=one
             end if
             ABI_DEALLOCATE(grhoe2)
             ABI_DEALLOCATE(grhocore2)
!            Compute contribution to annihilation rates
             ABI_ALLOCATE(ff,(mesh_size))
             if (option/=1) rhoarr1(:)=rhoarr1_j(:)
             do ii=1,4
               if (ii==1) ff(1:mesh_size)=rhoarr1(1:mesh_size)*rhoarr1_ep(1:mesh_size) &
&               *gamma(1:mesh_size,1)*pawrad(itypat)%rad(1:mesh_size)**2
               if (ii==2) ff(1:mesh_size)=rhoarr1(1:mesh_size)*rhoarr1_ep(1:mesh_size) &
&               *gamma(1:mesh_size,2)*pawrad(itypat)%rad(1:mesh_size)**2
               if (electronpositron%particle==EP_ELECTRON) then
                 if (ii==3) ff(1:mesh_size)=rhoarr1(1:mesh_size)*rhocore(1:mesh_size) &
&                 *gamma(1:mesh_size,1)*pawrad(itypat)%rad(1:mesh_size)**2
                 if (ii==4) ff(1:mesh_size)=rhoarr1(1:mesh_size)*rhocore(1:mesh_size) &
&                 *pawrad(itypat)%rad(1:mesh_size)**2
               else
                 if (ii==3) ff(1:mesh_size)=rhoarr1_ep(1:mesh_size)*rhocore(1:mesh_size) &
&                 *gamma(1:mesh_size,1)*pawrad(itypat)%rad(1:mesh_size)**2
                 if (ii==4) ff(1:mesh_size)=rhoarr1_ep(1:mesh_size)*rhocore(1:mesh_size) &
&                 *pawrad(itypat)%rad(1:mesh_size)**2
               end if
               call simp_gen(intg,ff,pawrad(itypat))
               intg=intg*pawang%angwgth(ipt)*four_pi
               if (ii==1) lambda_paw         =lambda_paw         +lsign(iloop)*intg
               if (ii==2) lambda_paw_ipm     =lambda_paw_ipm     +lsign(iloop)*intg
               if (ii==3) lambda_core_paw    =lambda_core_paw    +lsign(iloop)*intg
               if (ii==4) lambda_core_paw_ipm=lambda_core_paw_ipm+lsign(iloop)*intg
             end do
             ABI_DEALLOCATE(ff)
           end do ! ipt
           ABI_DEALLOCATE(gamma)
           ABI_DEALLOCATE(rhoarr1)
           ABI_DEALLOCATE(rhoarr1_ep)
           if (option/=1) then
             ABI_DEALLOCATE(rhoarr1_j)
           end if

!          Second formalism: use (l,m) moments for densities
         else if (dtset%pawxcdev/=0) then

!          Build densities
           ABI_ALLOCATE(gammam,(mesh_size,2,lm_size))
           ABI_ALLOCATE(rhotot,(mesh_size,lm_size))
           ABI_ALLOCATE(rhotot_ep,(mesh_size,lm_size))
           ABI_ALLOCATE(rhosph,(mesh_size))
           ABI_ALLOCATE(rhosph_ep,(mesh_size))
           if (option/=1) then
             ABI_ALLOCATE(rhotot_j,(mesh_size,lm_size))
             ABI_ALLOCATE(rhosph_j,(mesh_size))
           end if
           if (usecore==0) rhocore(:)=zero
           if (iloop==1) then
             rhotot   (:,:)=rho1   (:,:,1)
             rhotot_ep(:,:)=rho1_ep(:,:,1)
             if (option/=1) rhotot_j (:,:)=rho1_j (:,:,1)
             if (usecore==1) rhocore(:)=pawtab(itypat)%coredens(:)
           else
             if (include_nhat_in_gamma) then
               rhotot   (:,:)=trho1   (:,:,1)+nhat1   (:,:,1)
               rhotot_ep(:,:)=trho1_ep(:,:,1)+nhat1_ep(:,:,1)
               if (option/=1) rhotot_j (:,:)=trho1_j (:,:,1)+nhat1_j (:,:,1)
             else
               rhotot   (:,:)=trho1   (:,:,1)
               rhotot_ep(:,:)=trho1_ep(:,:,1)
               if (option/=1) rhotot_j (:,:)=trho1_j (:,:,1)
             end if
             if (usecore==1) rhocore(:)=pawtab(itypat)%tcoredens(:,1)
           end if
           rhosph   (:)=rhotot   (:,1)/sqfpi
           rhosph_ep(:)=rhotot_ep(:,1)/sqfpi
           if (option/=1) rhosph_j (:)=rhotot_j(:,1)/sqfpi
!          Make spherical densities positive
           if (electronpositron%particle==EP_ELECTRON) then
             call mkdenpos(iwarnp,mesh_size,1,1,rhosph   ,dtset%xc_denpos)
             call mkdenpos(iwarn ,mesh_size,1,1,rhosph_ep,dtset%xc_denpos)
           else if (electronpositron%particle==EP_POSITRON) then
             call mkdenpos(iwarn ,mesh_size,1,1,rhosph   ,dtset%xc_denpos)
             call mkdenpos(iwarnp,mesh_size,1,1,rhosph_ep,dtset%xc_denpos)
             if (option/=1) then
               call mkdenpos(iwarnp,mesh_size,1,1,rhosph_j,dtset%xc_denpos)
             end if
           end if
!          Need gradients of electronic densities for GGA
           ABI_ALLOCATE(grhoe2,(ngr))
           ABI_ALLOCATE(grhocore2,(ngr))
           if (ngr>0) then
             if (electronpositron%particle==EP_ELECTRON) then
               call nderiv_gen(grhoe2,rhosph_ep,pawrad(itypat))
             else if (electronpositron%particle==EP_POSITRON) then
               call nderiv_gen(grhoe2,rhosph,pawrad(itypat))
             end if
             grhoe2(:)=grhoe2(:)**2
             if (usecore==1) then
               call nderiv_gen(grhocore2,rhocore,pawrad(itypat))
               grhocore2(:)=grhocore2(:)**2
             end if
           end if
!          Compute Gamma for (rho-,rho+),
!          (rho- +drho-,rho+), (rho- -drho-,rho+),
!          (rho-,rho+ +drho+), (rho-,rho+ -drho+),
!          (rho- +drho-,rho+ +drho+), (rho- -drho-,rho+ -drho+)
!          Do a seven steps loop
           ABI_ALLOCATE(gam_,(mesh_size,2,7))
           ABI_ALLOCATE(rho_,(mesh_size))
           ABI_ALLOCATE(rho_ep_,(mesh_size))
           ABI_ALLOCATE(rhocor_,(mesh_size))
           ABI_ALLOCATE(grho2_,(ngr))
           ABI_ALLOCATE(grhocor2_,(ngr))
           do ii=1,7
!            Apply delta to get perturbed densities
             rho_(:)=rhosph(:);rho_ep_(:)=rhosph_ep(:);if (usecore==1) rhocor_(:)=rhocore(:)
             if (ngr>0) grho2_(:)=grhoe2(:)
             if (ngr>0) grhocor2_(:)=grhocore2(:)
             if (ii==2.or.ii==4.or.ii==6) fact=(one+delta)
             if (ii==3.or.ii==5.or.ii==7) fact=(one-delta)
             fact2=fact**2
             if (ii==2.or.ii==3.or.ii==6.or.ii==7) then
               rho_(:)=fact*rho_(:)
               if (electronpositron%particle==EP_POSITRON) then
                 if (ngr>0) grho2_(:)=fact2*grho2_(:)
                 if (usecore==1)rhocor_(:)=fact*rhocor_(:)
                 if (ngr>0.and.usecore==1) grhocor2_(:)=fact2*grhocor2_(:)
               end if
             end if
             if (ii==4.or.ii==5.or.ii==6.or.ii==7) then
               rho_ep_(:)=fact*rho_ep_(:)
               if (electronpositron%particle==EP_ELECTRON) then
                 if (ngr>0) grho2_(:)=fact2*grho2_(:)
                 if (usecore==1)rhocor_(:)=fact*rhocor_(:)
                 if (ngr>0.and.usecore==1) grhocor2_(:)=fact2*grhocor2_(:)
               end if
             end if
!            Compute gamma for these perturbed densities
             if (option==1.or.option==2) then
               if (electronpositron%particle==EP_ELECTRON) then
                 call gammapositron(gam_(:,:,ii),grhocor2_,grho2_,igamma(igam),ngr,mesh_size,rhocor_,rho_ep_,rho_,usecore)
               else if (electronpositron%particle==EP_POSITRON) then
                 call gammapositron(gam_(:,:,ii),grhocor2_,grho2_,igamma(igam),ngr,mesh_size,rhocor_,rho_,rho_ep_,usecore)
               end if
             else
               gam_(:,:,:)=one
             end if
           end do ! end loop ii=1,7

           ABI_DEALLOCATE(rhocor_)
           ABI_DEALLOCATE(grho2_)
           ABI_DEALLOCATE(grhocor2_)
           ABI_DEALLOCATE(grhoe2)
           ABI_DEALLOCATE(grhocore2)
           rho_   (:)=rhosph   (:);if (electronpositron%particle==EP_POSITRON.and.usecore==1) rho_   (:)=rho_   (:)+rhocore(:)
           rho_ep_(:)=rhosph_ep(:);if (electronpositron%particle==EP_ELECTRON.and.usecore==1) rho_ep_(:)=rho_ep_(:)+rhocore(:)
!          Compute numerical first and second derivatives of Gamma
!          d1gam(1) = dgam/drho+ (particle=ELECTRON), dgam/drho- (particle=POSITRON)
!          d1gam(2) = dgam/drho- (particle=ELECTRON), dgam/drho+ (particle=POSITRON)
           ABI_ALLOCATE(d1gam,(mesh_size,2,2))
           d1gam(:,:,:)=zero
           do ir=1,mesh_size
             if (rho_     (ir)>tol14) d1gam(ir,1,1)=(gam_(ir,1,2)-gam_(ir,1,3))*half/(delta*rho_     (ir))
             if (rhosph   (ir)>tol14) d1gam(ir,2,1)=(gam_(ir,2,2)-gam_(ir,2,3))*half/(delta*rhosph   (ir))
             if (rho_ep_  (ir)>tol14) d1gam(ir,1,2)=(gam_(ir,1,4)-gam_(ir,1,5))*half/(delta*rho_ep_  (ir))
             if (rhosph_ep(ir)>tol14) d1gam(ir,2,2)=(gam_(ir,2,4)-gam_(ir,2,5))*half/(delta*rhosph_ep(ir))
           end do

!          d2gam(1) = d2gam/drho+_drho+ (particle=ELECTRON), dgam/drho-_drho- (particle=POSITRON)
!          d2gam(2) = d2gam/drho-_drho+ (particle=ELECTRON), dgam/drho+_drho- (particle=POSITRON)
!          d2gam(3) = d2gam/drho-_drho- (particle=ELECTRON), dgam/drho+_drho+ (particle=POSITRON)
           ABI_ALLOCATE(d2gam,(mesh_size,2,3))
           d2gam(:,:,:)=zero
           do ir=1,mesh_size
             if (rho_  (ir)>tol14) d2gam(ir,1,1)=(gam_(ir,1,2)+gam_(ir,1,3)-two*gam_(ir,1,1))/(delta*rho_  (ir))**2
             if (rhosph(ir)>tol14) d2gam(ir,2,1)=(gam_(ir,2,2)+gam_(ir,2,3)-two*gam_(ir,2,1))/(delta*rhosph(ir))**2
             if (rho_ep_(ir)>tol14) then
               d2gam(ir,1,3)=(gam_(ir,1,4)+gam_(ir,1,5)-two*gam_(ir,1,1))/(delta*rho_ep_(ir))**2
               if (rho_(ir)>tol14) then
                 d2gam(ir,1,2)=(gam_(ir,1,6)+gam_(ir,1,7)+two*gam_(ir,1,1) &
&                 -gam_(ir,1,2)-gam_(ir,1,3)-gam_(ir,1,4)-gam_(ir,1,5)) &
&                 *half/(delta*rho_(ir))/(delta*rho_ep_(ir))
               end if
             end if
             if (rhosph_ep(ir)>tol14) then
               d2gam(ir,2,3)=(gam_(ir,2,4)+gam_(ir,2,5)-two*gam_(ir,2,1))/(delta*rhosph_ep(ir))**2
               if (rhosph(ir)>tol14) then
                 d2gam(ir,2,2)=(gam_(ir,2,6)+gam_(ir,2,7)+two*gam_(ir,2,1) &
&                 -gam_(ir,2,2)-gam_(ir,2,3)-gam_(ir,2,4)-gam_(ir,2,5)) &
&                 *half/(delta*rhosph(ir))/(delta*rhosph_ep(ir))
               end if
             end if
           end do
           ABI_DEALLOCATE(rho_)
           ABI_DEALLOCATE(rho_ep_)
!          Compute useful sums of densities
           ABI_ALLOCATE(v1sum,(mesh_size,3))
           if ( dtset%pawxcdev>=2)  then
             ABI_ALLOCATE(v2sum,(mesh_size,lm_size,3))
           else
             ABI_ALLOCATE(v2sum,(0,0,0))
           end if
           rhotot(:,1)=sqfpi*rhosph(:);rhotot_ep(:,1)=sqfpi*rhosph_ep(:)
           call pawxcsum(1,1,1,lmselect,lmselect_ep,lm_size,mesh_size,3,dtset%pawxcdev,&
&           pawang,rhotot,rhotot_ep,v1sum,v2sum)
!          Compute final developpment of gamma moments
           gammam(:,:,:)=zero
           gammam(:,:,1)=gam_(:,:,1)*sqfpi
           gammam(:,1,1)=gammam(:,1,1)+(d2gam(:,1,2)*v1sum(:,2) &
&           +half*(d2gam(:,1,1)*v1sum(:,1)+d2gam(:,1,3)*v1sum(:,3)))/sqfpi
           gammam(:,2,1)=gammam(:,2,1)+(d2gam(:,2,2)*v1sum(:,2) &
&           +half*(d2gam(:,2,1)*v1sum(:,1)+d2gam(:,2,3)*v1sum(:,3)))/sqfpi
           do ilm=2,lm_size
             if (lmselect(ilm)) then
               gammam(:,1,ilm)=gammam(:,1,ilm)+d1gam(:,1,1)*rhotot(:,ilm)
               gammam(:,2,ilm)=gammam(:,2,ilm)+d1gam(:,2,1)*rhotot(:,ilm)
             end if
             if (lmselect_ep(ilm)) then
               gammam(:,1,ilm)=gammam(:,1,ilm)+d1gam(:,1,2)*rhotot_ep(:,ilm)
               gammam(:,2,ilm)=gammam(:,2,ilm)+d1gam(:,2,2)*rhotot_ep(:,ilm)
             end if
           end do
           if (dtset%pawxcdev>1) then
             do ilm=2,lm_size
               gammam(:,1,ilm)=gammam(:,1,ilm)+d2gam(:,1,2)*v2sum(:,ilm,2) &
&               +half*(d2gam(:,1,1)*v2sum(:,ilm,1)+d2gam(:,1,3)*v2sum(:,ilm,3))
               gammam(:,2,ilm)=gammam(:,2,ilm)+d2gam(:,2,2)*v2sum(:,ilm,2) &
&               +half*(d2gam(:,2,1)*v2sum(:,ilm,1)+d2gam(:,2,3)*v2sum(:,ilm,3))
             end do
           end if
           ABI_DEALLOCATE(gam_)
           ABI_DEALLOCATE(d1gam)
           ABI_DEALLOCATE(d2gam)
           ABI_DEALLOCATE(v1sum)
           ABI_DEALLOCATE(v2sum)

!          Compute contribution to annihilation rate

!          In state dependent scheme for Doppler, replace electronic density with one state
           if (option/=1) then
             rhotot(:,:) = rhotot_j(:,:)
             rhosph  (:) = rhosph_j  (:)
           end if

           ABI_ALLOCATE(gg,(mesh_size,4))
           gg=zero
           ABI_ALLOCATE(rhoarr1,(mesh_size))
           ABI_ALLOCATE(rhoarr2,(mesh_size))
           do ilm=1,lm_size
             do ilm1=1,lm_size
               if (lmselect(ilm1)) then
                 if (ilm1==1) rhoarr1(:)=sqfpi*rhosph(:)
                 if (ilm1/=1) rhoarr1(:)=rhotot(:,ilm1)
                 do ilm2=1,lm_size
                   if (lmselect_ep(ilm2)) then
                     if (ilm2==1) rhoarr2(:)=sqfpi*rhosph_ep(:)
                     if (ilm2/=1) rhoarr2(:)=rhotot_ep(:,ilm2)
                     if (ilm1>=ilm2) then
                       isel=pawang%gntselect(ilm,ilm2+ilm1*(ilm1-1)/2)
                     else
                       isel=pawang%gntselect(ilm,ilm1+ilm2*(ilm2-1)/2)
                     end if
                     if (isel>0) then
                       fact=pawang%realgnt(isel)
                       gg(:,1)=gg(:,1)+fact*rhoarr1(:)*rhoarr2(:)*gammam(:,1,ilm)
                       gg(:,2)=gg(:,2)+fact*rhoarr1(:)*rhoarr2(:)*gammam(:,2,ilm)
                     end if
                   end if
                 end do
               end if
             end do
           end do
           ABI_DEALLOCATE(rhoarr1)
           ABI_DEALLOCATE(rhoarr2)
           if (electronpositron%particle==EP_ELECTRON) then
             do ilm=1,lm_size
               if (lmselect(ilm)) gg(:,3)=gg(:,3)+rhotot(:,ilm)*rhocore(:)*gammam(:,1,ilm)
             end do
             gg(:,4)=sqfpi*rhotot(:,1)*rhocore(:)
           else if (electronpositron%particle==EP_POSITRON) then
             do ilm=1,lm_size
               if (lmselect_ep(ilm)) gg(:,3)=gg(:,3)+rhotot_ep(:,ilm)*rhocore(:)*gammam(:,1,ilm)
             end do
             gg(:,4)=sqfpi*rhotot_ep(:,1)*rhocore(:)
           end if
           do ii=1,4
             gg(1:mesh_size,ii)=gg(1:mesh_size,ii)*pawrad(itypat)%rad(1:mesh_size)**2
             call simp_gen(intg,gg(:,ii),pawrad(itypat))
             if (ii==1) lambda_paw         =lambda_paw         +lsign(iloop)*intg
             if (ii==2) lambda_paw_ipm     =lambda_paw_ipm     +lsign(iloop)*intg
             if (ii==3) lambda_core_paw    =lambda_core_paw    +lsign(iloop)*intg
             if (ii==4) lambda_core_paw_ipm=lambda_core_paw_ipm+lsign(iloop)*intg
           end do
           ABI_DEALLOCATE(gg)
           ABI_DEALLOCATE(gammam)
           ABI_DEALLOCATE(rhotot)
           ABI_DEALLOCATE(rhotot_ep)
           ABI_DEALLOCATE(rhosph)
           ABI_DEALLOCATE(rhosph_ep)
           if (option/=1) then
             ABI_DEALLOCATE(rhotot_j)
             ABI_DEALLOCATE(rhosph_j)
           end if

         end if ! dtset%pawxcdev

         ABI_DEALLOCATE(rhocore)

       end do ! iloop

       ABI_DEALLOCATE(rho1)
       ABI_DEALLOCATE(trho1)
       ABI_DEALLOCATE(rho1_ep)
       ABI_DEALLOCATE(trho1_ep)
       ABI_DEALLOCATE(rho1_j)
       ABI_DEALLOCATE(trho1_j)
       ABI_DEALLOCATE(nhat1)
       ABI_DEALLOCATE(nhat1_ep)
       ABI_DEALLOCATE(nhat1_j)
       ABI_DEALLOCATE(lmselect)
       ABI_DEALLOCATE(lmselect_ep)
       ABI_DEALLOCATE(lmselect_dum)

     end do ! iatom

!    Reduction in case of distribution over atomic sites
     if (mpi_enreg%nproc_atom>1) then
       mpibuf(1)=lambda_paw     ;mpibuf(2)=lambda_paw_ipm
       mpibuf(3)=lambda_core_paw;mpibuf(4)=lambda_core_paw_ipm
       call xmpi_sum(mpibuf,mpi_enreg%comm_atom,ierr)
       lambda_paw=mpibuf(1)     ;lambda_paw_ipm=mpibuf(2)
       lambda_core_paw=mpibuf(3);lambda_core_paw_ipm=mpibuf(4)
     end if

!    Add plane-wave and PAW contributions to annihilation rates
     lambda         =lambda         +lambda_paw
     lambda_ipm     =lambda_ipm     +lambda_paw_ipm
     lambda_core    =lambda_core    +lambda_core_paw
     lambda_core_ipm=lambda_core_ipm+lambda_core_paw_ipm
   end if ! dtset%usepaw


!  Convert into proper units and print
!  ---------------------------------------------------------------------------------------

!  Sum valence and core contributions to annihilation rates
   lambda        =lambda        +lambda_core
   lambda_ipm    =lambda_ipm    +lambda_core_ipm
   if (dtset%usepaw==1) then
     lambda_paw    =lambda_paw    +lambda_core_paw
     lambda_paw_ipm=lambda_paw_ipm+lambda_core_paw_ipm
   end if

!  Set annihilation rate in proper unit (picosec.)
   units=pi*(one/InvFineStruct)**3/Time_Sec/1.e12_dp/electronpositron%posocc

   lambda         =lambda         *units
   lambda_ipm     =lambda_ipm     *units
   lambda_core    =lambda_core    *units
   lambda_core_ipm=lambda_core_ipm*units
   lifetime    =one/lambda
   lifetime_ipm=one/lambda_ipm
   electronpositron%lambda=lambda
   electronpositron%lifetime=lifetime
   if (dtset%usepaw==1) then
     lambda_paw         =lambda_paw         *units
     lambda_paw_ipm     =lambda_paw_ipm     *units
     lambda_core_paw    =lambda_core_paw    *units
     lambda_core_paw_ipm=lambda_core_paw_ipm*units
   end if
   rate=lambda-lambda_core-lambda_paw+lambda_core_paw
   rate_paw=lambda_paw-lambda_core_paw
!  Print life time and additional information
   if (option==1) then
     if (igam==1) then
       write(msg,'(a,80("-"),2a)') ch10,ch10,' Results for electron-positron annihilation:'
       call wrtout(ab_out,msg,'COLL')
       call wrtout(std_out,msg,'COLL')
     end if
     if (ngamma>1.and.igam==1) then
       write(msg,'(a,i1,a)') ch10,ngamma,&
&       ' computations of positron lifetime have been performed (with different enhancement factors).'
       call wrtout(ab_out,msg,'COLL')
       call wrtout(std_out,msg,'COLL')
     end if
     if (ngamma>1) then
       write(msg,'(2a,i1)') ch10,"########## Lifetime computation ",igam
       call wrtout(ab_out,msg,'COLL')
       call wrtout(std_out,msg,'COLL')
     end if
     if (abs(electronpositron%ixcpositron)==1) then
       write(msg,'(4a)') ch10,' # Zero-positron density limit of Arponen and Pajanne provided by Boronski & Nieminen',&
&       ch10,'   Ref.: Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986)' ! [[cite:Boronski1986]]
     else if (electronpositron%ixcpositron==11) then
       write(msg,'(4a)') ch10,' # Zero-positron density limit of Arponen and Pajanne fitted by Sterne & Kaiser',&
&       ch10,'   Ref.: P.A. Sterne and J.H. Kaiser, Phys. Rev. B 43, 13892 (1991)' ! [[cite:Sterne1991]]
     else if (electronpositron%ixcpositron==2) then
       write(msg,'(4a)') ch10,' # Electron-positron correlation provided by Puska, Seitsonen, and Nieminen',&
&       ch10,'   Ref: M.J. Puska, A.P. Seitsonen and R.M. Nieminen, Phys. Rev. B 52, 10947 (1994)' !  [[cite:Puska1994]]
     else if (electronpositron%ixcpositron==3) then
       write(msg,'(8a)') ch10,' # Zero-positron density limit of Arponen and Pajanne provided by Boronski & Nieminen',&
&       ch10,'   + GGA corrections',&
&       ch10,'   Ref.: Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986)',& ! [[cite:Boronski1986]]
&       ch10,'         B. Barbiellini, M.J. Puska, T. Torsti and R.M.Nieminen, Phys. Rev. B 51, 7341 (1995)' ! [[cite:Barbiellini1995]]
     else if (electronpositron%ixcpositron==31) then
       write(msg,'(8a)') ch10,' # Zero-positron density limit of Arponen and Pajanne fitted by Sterne & Kaiser',&
&       ch10,'   + GGA corrections',&
&       ch10,'   Ref.: P.A. Sterne and J.H. Kaiser, Phys. Rev. B 43, 13892 (1991)',& ! [[cite:Sterne1991]]
&       ch10,'         B. Barbiellini, M.J. Puska, T. Torsti and R.M. Nieminen, Phys. Rev. B 51, 7341 (1995)' ! [[cite:Barbiellini1995]]
     end if
     call wrtout(ab_out,msg,'COLL')
     call wrtout(std_out,  msg,'COLL')
     if (igamma(igam)==0) then
       write(msg,'(a)')       ' # Enhancement factor set to one (test)'
     else if (igamma(igam)==1) then
       write(msg,'(3a)')      ' # Enhancement factor of Boronski & Nieminen',&
&       ch10,'   Ref.: Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986)' ! [[cite:Boronski1986]]
     else if (igamma(igam)==2) then
       write(msg,'(3a)')      ' # Enhancement factor of Boronski & Nieminen IN THE RPA LIMIT',&
&       ch10,'   Ref.: Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986)' ! [[cite:Boronski1986]]
     else if (igamma(igam)==3) then
       write(msg,'(3a)')      ' # Enhancement factor of Sterne & Kaiser',&
&       ch10,'   Ref.: P.A. Sterne and J.H. Kaiser, Phys. Rev. B 43, 13892 (1991)' ! [[cite:Sterne1991]]
     else if (igamma(igam)==4) then
       write(msg,'(3a)')      ' # Enhancement factor of Puska, Seitsonen, and Nieminen',&
&       ch10,'   Ref.: M.J. Puska, A.P. Seitsonen and R.M. Nieminen, Phys. Rev. B 52, 10947 (1994)' !  [[cite:Puska1994]]
     end if
     call wrtout(ab_out,msg,'COLL')
     call wrtout(std_out,msg,'COLL')
     write(msg, '(4(2a,es16.8))' ) ch10,&
&     ' Positron lifetime                         (ps)   =',lifetime    ,ch10,&
&     ' Positron lifetime with IPM for core elec. (ps)   =',lifetime_ipm,ch10,&
&     ' Annihilation rate                         (ns-1) =',lambda    *1000._dp,ch10,&
&     ' Annihilation rate with IPM for core elec. (ns-1) =',lambda_ipm*1000._dp
     call wrtout(ab_out,msg,'COLL')
     call wrtout(std_out,msg,'COLL')
     write(msg,'(2a,5(2a,es16.8))' ) ch10,&
&     ' Annihilation rate core/valence decomposition:',ch10,&
&     '   Core    contribution to ann.rate          (ns-1) =', lambda_core                 *1000._dp,ch10,&
&     '   Valence contribution to ann.rate          (ns-1) =',(lambda-lambda_core)         *1000._dp,ch10,&
&     '   Core    contribution to ann.rate with IPM (ns-1) =', lambda_core_ipm             *1000._dp,ch10,&
&     '   Valence contribution to ann.rate with IPM (ns-1) =',(lambda_ipm-lambda_core_ipm) *1000._dp
     call wrtout(ab_out,msg,'COLL')
     call wrtout(std_out,msg,'COLL')
     if (dtset%usepaw==1) then
       write(msg, '(2a,6(2a,es16.8))' ) ch10,&
&       ' Annihilation rate PAW decomposition:',ch10,&
&       '   Plane-wave contribution to ann.rate          (ns-1) =',(lambda-lambda_paw)*1000._dp,ch10,&
&       '   Plane-wave valence contribution to ann.rate  (ns-1) =',(lambda-lambda_paw-lambda_core+lambda_core_paw)*1000._dp,ch10,&
&       '   On-site core contribution to ann.rate        (ns-1) =', lambda_core_paw*1000._dp,ch10,&
&       '   On-site valence contribution to ann.rate     (ns-1) =',(lambda_paw-lambda_core_paw)*1000._dp,ch10,&
&       '   Plane-wave contribution to ann.rate with IPM (ns-1) =',(lambda_ipm-lambda_paw_ipm)*1000._dp,ch10,&
&       '   Plane-wave core contrb. to ann.rate with IPM (ns-1) =',(lambda_core_ipm-lambda_core_paw_ipm)*1000._dp
       call wrtout(ab_out,msg,'COLL')
       call wrtout(std_out,msg,'COLL')
     end if
     if (dtset%usepaw==0.and.igam==ngamma) then ! These tests are not relevant with PAW
       write(msg, '(2a,3(2a,es16.8))' ) ch10,&
&       ' ########## Some checks, for testing purpose:',ch10,&
&       '   Number of core electrons      =',nbec,ch10,&
&       '   Number of valence electrons   =',nbev,ch10,&
&       '   Number of positrons           =',nbp
       call wrtout(ab_out,msg,'COLL')
       call wrtout(std_out,msg,'COLL')
     end if
   end if !end if option
 end do ! Big loop on igam

 if (option==1) then
   write(msg, '(3a)' ) ch10,'      (*) IPM=Independent particle Model',ch10
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
 end if !end if option

!Deallocate memory
 ABI_DEALLOCATE(igamma)
 if (dtset%usepaw==1.and.(.not.include_nhat_in_gamma)) then
   ABI_DEALLOCATE(rhor_)
   ABI_DEALLOCATE(rhor_ep_)
 end if

 DBG_EXIT("COLL")

end subroutine poslifetime
!!***

!!****f* ABINIT/posdoppler
!! NAME
!! posdoppler
!!
!! FUNCTION
!! Calculate the momentum distribution annihilating electrons-positron (Doppler broadening)
!!
!! INPUTS
!!  cg(2,mcg)=planewave coefficients of wavefunctions.
!!  cprj(natom,mcprj)= <p_lmn|Cnk> coefficients for each WF |Cnk>
!!                     and each |p_lmn> non-local projector
!!  Cryst<Crystal_structure> = Info on unit cell and its symmetries
!!  dimcprj(natom)=array of dimensions of array cprj (not ordered)
!!  dtfil <type(datafiles_type)>=variables related to files
!!   | unpaw=unit number for temporary PAW files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | istwfk=input option=1 parameter that describes the storage of wfs
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs for the "coarse" grid
!!   | mkmem=number of k points treated by this node.
!!   | mpw=maximum dimensioned size of npw
!!   | natom=number of atoms
!!   | nband=number of bands at each k point
!!   | ngfft=contain all needed information about 3D FFT (coarse grid)
!!   | nkpt=number of k points
!!   | nspden=number of spin-density components
!!   | nspinor=number of spinorial components of the wavefunctions
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | usepaw=flag for PAW
!!   | use_gpu_cuda=flag for Cuda use
!!   | wtk(=weights associated with various k points
!!  filpsp(ntypat)=name(s) of the pseudopotential file(s)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mpi_enreg= informations about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  n3xccc= dimension of the xccc3d array (0 or nfft).
!!  nfft= number of FFT grid points
!!  ngfft(18)= contain all needed information about 3D FFT
!!  nhat(nfft,nspden)=charge compensation density (content depends on electronpositron%particle)
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  occ(mband*nkpt*nsppol)=occupancies for each band and k point
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(my_natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  rhor(nfft,nspden)=total electron/positron density (content depends on electronpositron%particle)
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!
!! TODO
!!  print a warning if the core wave function is not localized in the PAW sphere
!!  implement PAW on-site contribution for state-independent scheme
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      bandfft_kpt_destroy,bandfft_kpt_mpi_recv,bandfft_kpt_mpi_send
!!      bandfft_kpt_reset,destroy_mpi_enreg,fourdp,fourwf,gammapositron_fft
!!      initmpi_seq,initylmr,mkdenpos,pawaccrhoij,pawcprj_alloc,pawcprj_bcast
!!      pawcprj_copy,pawcprj_free,pawcprj_get,pawcprj_mpi_recv,pawcprj_mpi_send
!!      pawpsp_read_corewf,pawrhoij_alloc,pawrhoij_free,pawrhoij_gather
!!      pawrhoij_nullify,poslifetime,posratecore,prep_fourwf,ptabs_fourdp,sbf8
!!      set_mpi_enreg_fft,simp_gen,sphereboundary,pawrhoij_symrhoij,unset_mpi_enreg_fft
!!      wffclose,wffopen,wrtout,xderivewrecend,xderivewrecinit,xderivewrite
!!      xmoveoff,xmpi_bcast,xmpi_recv,xmpi_send,xmpi_sum
!!
!! SOURCE

!Macro to go from row-column indexing to combined indexing
#define RCC(glmn,hlmn) max(glmn,hlmn)*(max(glmn,hlmn)-1)/2+min(glmn,hlmn)
!Macro to go from l,m angular momentum indexing to combined indexing
#define LMC(lval,mval) lval*lval+lval+mval+1

subroutine posdoppler(cg,cprj,Crystal,dimcprj,dtfil,dtset,electronpositron,&
&                     filpsp,kg,mcg,mcprj,mpi_enreg,my_natom,&
&                     n3xccc,nfft,ngfft,nhat,npwarr,occ,pawang,pawrad,&
&                     pawrhoij,pawtab,rhor,xccc3d)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mcg,mcprj,my_natom,n3xccc,nfft
 type(crystal_t) :: Crystal
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer :: electronpositron
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: dimcprj(dtset%natom)
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),ngfft(18),npwarr(dtset%nkpt)
 real(dp),intent(in) :: nhat(nfft,dtset%nspden*dtset%usepaw),xccc3d(n3xccc)
 real(dp),intent(in),target :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in),target :: rhor(nfft,dtset%nspden)
 real(dp),intent(inout),target :: cg(2,mcg)
 character(len=fnlen),intent(in) :: filpsp(dtset%ntypat)
 type(pawcprj_type),target :: cprj(dtset%natom,mcprj)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*dtset%usepaw)
 type(pawrhoij_type),intent(in),target :: pawrhoij(my_natom*dtset%usepaw)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer :: accessfil,basis_size,bandpp,bdtot_index,bdtot_index_pos,blocksize,cplex,cplex_rhoij
 integer :: glmij,i0lmn,i1,i2,i3,iat,iatm,iatom
 integer :: ib,ib_cprj,ib_cprj_pos,ib_pos,ibg,ibg_pos
 integer :: iblock,iblock_pos,ibpp,ibpp_pos
 integer :: icg,icg_pos,id1,id2,id3,ierr,ig1,ig2,ig3,igamma,ii,ikg,ikg_pos,ikpt
 integer :: ikpt_pos,il,ilm,ilmn,iln,indx,indx0,iorder_cprj,iproc,ir,isppol,isppol_pos,istwf_k
 integer :: istwf_k_pos,itypat,iwarn,iwavef,iwavef_pos,j2,j3,jj,jkpt,jl,jlm,jlmn,jln
 integer :: klm,kln,klmn,l_size,ll,llmax,llmin,lm,lmn_size,lmn2_size
 integer :: mband_cprj,mband_cprj_pos,mcg_pos
 integer :: mcprj_k,mcprj_k_pos,me_band,me_fft,me_kpt,me_kptband
 integer :: mesh_size,meshsz,mm,my_ngrid,my_nspinor,my_nsppol,my_n2,n1,n2,n3,n4,n5,n6
 integer :: nband_cprj_eff_pos,nband_cprj_k,nband_cprj_k_pos
 integer :: nband_eff_pos,nband_k,nband_k_pos
 integer :: nblock_band,nblock_band_eff_pos,nkpt
 integer :: nproc_band,nproc_fft,nproc_kpt,nproc_kptband,npw_k,npw_k_pos
 integer :: nspden_rhoij,option,tag,unit_doppler
 integer :: tim_fourdp=0,tim_fourwf=-36
 integer :: ylmr_normchoice,ylmr_npts,ylmr_option
 logical,parameter :: include_nhat_in_gamma=.false.,state_dependent=.true.
 logical,parameter :: kgamma_only_positron=.true.,wf_conjugate=.false.
 logical :: cprj_paral_band,ex,mykpt,mykpt_pos,usetimerev
 real(dp) :: arg,bessarg,cpi,cpr,cp11,cp12,cp21,cp22,gammastate,intg
 real(dp) :: lambda_v1,lambda_v2,lambda_core,lambda_pw,occ_el,occ_pos
 real(dp) :: pnorm,pr,rate,rate_ipm,ratec,ratec_ipm,rate_paw,rate_paw_ipm
 real(dp) :: scale_,units_,weight,weight_pos,wf_fact,wtk_k,wtk_k_pos,vec
 character(len=fnlen) :: filename,filename_dop
 character(len=500) :: msg
 type(bandfft_kpt_type),pointer :: bandfft_kpt_el,bandfft_kpt_pos
 type(MPI_type) :: mpi_enreg_seq
 type(wffile_type) :: wff
!arrays
 integer,allocatable :: gbound(:,:),gbound_pos(:,:),kg_k(:,:),kg_k_pos(:,:)
 integer,allocatable :: lcor(:),lmncmax(:),my_ffttab(:),my_gridtab(:),ncor(:),nphicor(:)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 logical,allocatable :: have_intc(:,:,:),have_rad(:,:)
 real(dp) :: buf(4),contrib(2),cp(2),cp_pos(2),expipr(2),pbn(3),pcart(3)
 real(dp) :: radsumnfftc(2),ylmgr(1,1,0),ylmr_nrm(1)
 real(dp),allocatable :: cwaveg(:,:),cwaveg_pos(:,:),cwaver(:),cwaver_pos(:),cwaver_pos_block(:)
 real(dp),allocatable :: cg_k_pos(:,:),cwaveaug(:,:,:,:),cwaveaug_pos(:,:,:,:)
 real(dp),allocatable :: denpot_dum(:,:,:),energycor(:),ff(:),fofgout_dum(:,:)
 real(dp),allocatable :: gamma(:,:),intc(:,:,:),j_bessel(:,:),jbes(:),mpibuf(:,:)
 real(dp),allocatable :: occ_k(:),occ_k_pos(:),pcart_k(:,:)
 real(dp),allocatable :: radint1(:,:),radint2(:,:),radint3(:,:)
 real(dp),allocatable :: radsumnfft1(:,:),radsumnfft2(:,:),radsumnfft3(:,:)
 real(dp),allocatable :: rho_contrib(:),rho_contrib_g(:,:)
 real(dp),allocatable :: rho_contrib_paw1(:,:),rho_contrib_paw2(:,:),rho_contrib_paw3(:,:)
 real(dp),allocatable :: rho_moment_v1(:,:),rho_moment_v2(:,:)
 real(dp),allocatable :: rho_moment_core(:,:),rho_moment_k(:),rho_moment_k2(:)
 real(dp),allocatable :: rho_pw(:,:),rhor_dop_el(:)
 real(dp),allocatable :: rhocorej(:),rhoe(:,:),rhop(:,:),ylmp(:)
 real(dp),pointer :: cg_pos_ptr(:,:),cg_ptr(:,:),occ_ptr(:),occ_pos_ptr(:)
 real(dp),pointer :: rhor_(:,:),rhor_ep_(:,:)
 complex(dpc) :: ifac ! (-i)^L mod 4
 complex(dpc),dimension(0:3) :: ilfac(0:3)=(/(1.0,0.0),(0.0,-1.0),(-1.0,0.0),(0.0,1.0)/)
 type(coeff1_type),allocatable :: gammastate_c(:)
 type(coeffi2_type),allocatable :: indlmncor(:)
 type(coeff2_type),allocatable :: phicor(:)
 type(coeff6_type),allocatable :: radsum1(:),radsum2(:),radsum3(:)
 type(coeff7_type),allocatable :: radsumc(:)
 type(pawcprj_type),allocatable :: cprj_k(:,:),cprj_k_pos(:,:),cprj_pos(:,:)
 type(pawcprj_type),pointer :: cprj_pos_ptr(:,:),cprj_ptr(:,:)
 type(pawrhoij_type),allocatable :: pawrhoij_dop_el(:)
 type(pawrhoij_type),pointer :: pawrhoij_ptr(:),pawrhoij_all(:),pawrhoij_ep_all(:)

! *************************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 if (.not.associated(electronpositron)) then
   msg='electronpositron variable must be associated!'
   MSG_BUG(msg)
 end if
 if (allocated(mpi_enreg%proc_distrb)) then
   do isppol=1,dtset%nsppol
     do ikpt=1,dtset%nkpt
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
       if (any(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)/=mpi_enreg%proc_distrb(ikpt,1,isppol))) then
         msg='proc_distrib cannot be distributed over bands!'
         MSG_BUG(msg)
       end if
     end do
   end do
 end if
 if (dtset%nspinor==2) then
   msg='Doppler broadening not available for spinorial wave functions (nspinor=2)!'
   MSG_BUG(msg)
 end if
 if (mcprj==0) then
   msg='<p|Psi> (cprj) datastructure must be kept in memory (see pawusecp input keyword)!'
   MSG_BUG(msg)
 end if
 if (dtset%usepaw==0) then
   write(msg,'(5a)') 'Momentum distribution of annihilating electron-positron pairs',ch10,&
&   'in the Norm-conserving Pseudopotential formalism is incomplete!',ch10,&
&   'No core contribution is included.'
   MSG_WARNING(msg)
 end if
 if (any(dtset%nband(:)/=dtset%nband(1))) then
   write(msg,'(a)') 'Number of bands has to be the same for all k-points!'
   MSG_BUG(msg)
 end if
 if (dtset%usepaw==1) then
   if (size(pawrhoij)/=mpi_enreg%my_natom) then
     write(msg,'(a)') 'wrong size for pawrhoij! '
     MSG_BUG(msg)
   end if
 end if

!Various initializations
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 n4=ngfft(4) ; n5=ngfft(5) ; n6=ngfft(6)
 id1=n1/2+2 ; id2=n2/2+2 ; id3=n3/2+2
 iorder_cprj=0 ; cplex=2 ; iwarn=1
 wf_fact=one;if (wf_conjugate) wf_fact=-one
 nkpt=dtset%nkpt

!Manage kpt/spin parallelism
 ABI_ALLOCATE(my_gridtab,(nkpt))
 my_gridtab=0
 do ii=1,nkpt
   if (any(mpi_enreg%my_isppoltab(:)==1)) my_gridtab(ii)=mpi_enreg%my_kpttab(ii)
 end do
 my_ngrid=count(my_gridtab(:)/=0)
 my_nsppol=sum(mpi_enreg%my_isppoltab(:))
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

!Parallel settings
 if (mpi_enreg%paral_kgb/=0) then
   nproc_kpt=mpi_enreg%nproc_kpt
   nproc_band=mpi_enreg%nproc_band
   nproc_fft=mpi_enreg%nproc_fft
   nproc_kptband=xmpi_comm_size(mpi_enreg%comm_kptband)
   me_kpt=mpi_enreg%me_kpt
   me_band=mpi_enreg%me_band
   me_fft=mpi_enreg%me_fft
   me_kptband=xmpi_comm_rank(mpi_enreg%comm_kptband)
   bandpp=mpi_enreg%bandpp
   my_n2=n2/nproc_fft
   accessfil=IO_MODE_FORTRAN;if(nproc_fft>1)accessfil=IO_MODE_MPI
 else
   nproc_kpt=mpi_enreg%nproc_kpt
   nproc_band=1;nproc_fft=1
   nproc_kptband=nproc_kpt
   me_band=0;me_fft=0
   me_kpt=mpi_enreg%me_kpt
   me_kptband=me_kpt
   bandpp=1 ; my_n2=n2
   accessfil=IO_MODE_FORTRAN
 end if
 blocksize=nproc_band*bandpp
 nblock_band=dtset%nband(1)/blocksize

!Select density according to nhat choice0
 if (dtset%usepaw==0.or.include_nhat_in_gamma) then
   rhor_ => rhor
   rhor_ep_ => electronpositron%rhor_ep
 else
   ABI_ALLOCATE(rhor_,(nfft,dtset%nspden))
   ABI_ALLOCATE(rhor_ep_,(nfft,dtset%nspden))
   rhor_=rhor-nhat
   rhor_ep_=electronpositron%rhor_ep-electronpositron%nhat_ep
 end if

!Select type(s) of enhancement factor
 igamma=0
 if (electronpositron%ixcpositron==-1) igamma=0
 if (electronpositron%ixcpositron== 1) igamma=2
 if (electronpositron%ixcpositron== 2) igamma=4
 if (electronpositron%ixcpositron== 3) igamma=2
 if (electronpositron%ixcpositron==11) igamma=3
 if (electronpositron%ixcpositron==31) igamma=3

!Select electronic and positronic states
 if (electronpositron%particle==EP_ELECTRON) then !we should not be in this case
   cg_ptr => electronpositron%cg_ep
   cprj_ptr => electronpositron%cprj_ep
   occ_ptr => electronpositron%occ_ep
   pawrhoij_ptr => electronpositron%pawrhoij_ep
   cg_pos_ptr => cg
   cprj_pos_ptr => cprj
   occ_pos_ptr => occ
 end if
 if (electronpositron%particle==EP_POSITRON) then
   cg_ptr => cg
   cprj_ptr => cprj
   occ_ptr => occ
   pawrhoij_ptr => pawrhoij
   cg_pos_ptr => electronpositron%cg_ep
   cprj_pos_ptr => electronpositron%cprj_ep
   occ_pos_ptr => electronpositron%occ_ep
 end if

!Determine if cprj datastructures are distributed over bands
 mband_cprj=size(cprj_ptr,2)/(my_nspinor*dtset%mkmem*dtset%nsppol)
 mband_cprj_pos=size(cprj_pos_ptr,2)/(my_nspinor*dtset%mkmem*dtset%nsppol)
 cprj_paral_band=(mband_cprj<dtset%mband)

!Get the distrib associated with the fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!===============================================================================
!================ Calculate the PAW on-site constants ==========================

 if (dtset%usepaw==1) then

   ylmr_normchoice = 0 ! input to initylmr are normalized
   ylmr_npts = 1 ! only 1 point to compute in initylmr
   ylmr_nrm(1) = one ! weight of normed point for initylmr
   ylmr_option = 1 ! compute only ylm's in initylmr

  !Prepare radial integral for PAW correction for each atom type
   ABI_DATATYPE_ALLOCATE(radsum1,(dtset%ntypat))
   ABI_DATATYPE_ALLOCATE(radsum2,(dtset%ntypat))
   ABI_DATATYPE_ALLOCATE(radsum3,(dtset%ntypat))
   ABI_DATATYPE_ALLOCATE(radsumc,(dtset%ntypat))

   ABI_DATATYPE_ALLOCATE(indlmncor,(dtset%ntypat))
   ABI_DATATYPE_ALLOCATE(phicor,(dtset%ntypat))
   ABI_DATATYPE_ALLOCATE(gammastate_c,(dtset%natom))
   ABI_ALLOCATE(nphicor,(dtset%ntypat))
   ABI_ALLOCATE(lmncmax,(dtset%ntypat))

!  Reading of core wave functions
   if (mpi_enreg%me_cell==0) then
     do itypat=1,dtset%ntypat
       filename=trim(filpsp(itypat))//'.corewf'
       inquire(file=filename,exist=ex)
       if (.not.ex) then
         write(unit=filename,fmt='(a,i1)') 'corewf.abinit',itypat
         inquire(file=filename,exist=ex)
         if (.not.ex) then
           msg='Core wave-functions file is missing!'
           MSG_ERROR(msg)
         end if
       end if
       call pawpsp_read_corewf(energycor,indlmncor(itypat)%value,lcor,lmncmax(itypat),&
&       ncor,nphicor(itypat),pawrad(itypat),phicor(itypat)%value,&
&       filename=filename)
!      The following arrays are not used anymore
       ABI_DEALLOCATE(energycor)
       ABI_DEALLOCATE(lcor)
       ABI_DEALLOCATE(ncor)
     end do
   end if
   if (mpi_enreg%nproc_cell>1) then
     call xmpi_bcast(indlmncor,0,mpi_enreg%comm_cell,ierr)
     call xmpi_bcast(phicor,0,mpi_enreg%comm_cell,ierr)
     call xmpi_bcast(nphicor,0,mpi_enreg%comm_cell,ierr)
     call xmpi_bcast(lmncmax,0,mpi_enreg%comm_cell,ierr)
   end if

   do itypat=1,dtset%ntypat

     mesh_size = pawtab(itypat)%mesh_size
     l_size = pawtab(itypat)%l_size
     lmn_size = pawtab(itypat)%lmn_size
     lmn2_size = pawtab(itypat)%lmn2_size
     basis_size = pawtab(itypat)%basis_size

     ABI_ALLOCATE(j_bessel,(mesh_size,l_size))
     ABI_ALLOCATE(ylmp,(l_size*l_size))
     ABI_ALLOCATE(have_intc,(l_size,basis_size,nphicor(itypat)))
     ABI_ALLOCATE(have_rad,(l_size,pawtab(itypat)%ij_size))
     ABI_ALLOCATE(intc,(l_size,basis_size,nphicor(itypat)))
     ABI_ALLOCATE(radint1,(l_size,pawtab(itypat)%ij_size))
     ABI_ALLOCATE(radint2,(l_size,pawtab(itypat)%ij_size))
     ABI_ALLOCATE(radint3,(l_size,pawtab(itypat)%ij_size))

     ABI_ALLOCATE(radsumc(itypat)%value,(2,lmn_size,lmncmax(itypat),n1,my_n2,n3,my_ngrid))
     ABI_ALLOCATE(radsum1(itypat)%value,(2,lmn2_size,n1,my_n2,n3,my_ngrid))
     ABI_ALLOCATE(radsum2(itypat)%value,(2,lmn2_size,n1,my_n2,n3,my_ngrid))
     ABI_ALLOCATE(radsum3(itypat)%value,(2,lmn2_size,n1,my_n2,n3,my_ngrid))
     radsumc(itypat)%value=zero
     radsum1(itypat)%value=zero
     radsum2(itypat)%value=zero
     radsum3(itypat)%value=zero

     ABI_ALLOCATE(jbes,(l_size))
     ABI_ALLOCATE(ff,(mesh_size))
     meshsz=pawrad(itypat)%int_meshsz
     if (meshsz>mesh_size) ff(meshsz+1:mesh_size)=zero

     indx=0;jkpt=0
     do ikpt=1,nkpt
       if (my_gridtab(ikpt)==0) cycle
       jkpt=jkpt+1
       do i3=1,n3
         ig3=i3-(i3/id3)*n3-1
         do i2=1,n2
           if (me_fft/=fftn2_distrib(i2)) cycle
           j2=ffti2_local(i2)
           indx=n1*(my_n2*(i3-1)+(j2-1))
           ig2=i2-(i2/id2)*n2-1
           do i1=1,n1
             ig1=i1-(i1/id1)*n1-1
             indx=indx+1;if (mod(indx-1,nproc_band)/=me_band) cycle

             pcart(:)=Crystal%gprimd(:,1)*real(ig1+dtset%kpt(1,ikpt))+&
&             Crystal%gprimd(:,2)*real(ig2+dtset%kpt(2,ikpt))+&
&             Crystal%gprimd(:,3)*real(ig3+dtset%kpt(3,ikpt))
             pnorm=dsqrt(dot_product(pcart,pcart))

             if (pnorm < tol12) then
               pbn(:) = zero
               ylmp(:) = zero
               ylmp(1) = 1.d0/sqrt(four_pi)
             else
               pbn(:) = pcart(:)/pnorm ! unit vector
               call initylmr(l_size,ylmr_normchoice,ylmr_npts,ylmr_nrm,ylmr_option,pbn(:),ylmp(:),ylmgr)
             end if

             pnorm=two_pi*pnorm ! re-normed for call to bessel
             do ir = 1, mesh_size
               bessarg = pnorm*pawrad(itypat)%rad(ir)
               call sbf8(l_size,bessarg,jbes)
               j_bessel(ir,:)=jbes(:)
             end do

!            ===== Core part =====
!            Need intc=\int phi phi_core jl (pr) dr

             have_intc(:,:,:)=.FALSE. ; intc(:,:,:)=zero

             do jlmn = 1,lmncmax(itypat)
               jln = indlmncor(itypat)%value(5,jlmn)
               jlm = indlmncor(itypat)%value(4,jlmn)
               jl  = indlmncor(itypat)%value(1,jlmn)
               do ilmn = 1,lmn_size
                 iln = pawtab(itypat)%indlmn(5,ilmn)
                 ilm = pawtab(itypat)%indlmn(4,ilmn)
                 il  = pawtab(itypat)%indlmn(1,ilmn)

                 llmin = abs(il-jl)
                 llmax = il+jl
                 klm = RCC(ilm,jlm)
                 do ll=llmin,llmax,2
                   ifac=ilfac(mod(ll,4))

                   if (.not.have_intc(ll+1,iln,jln)) then
                     ff(1:mesh_size)=(pawtab(itypat)%phi(1:mesh_size,iln)*phicor(itypat)%value(1:mesh_size,jln))&
&                     *j_bessel(1:mesh_size,ll+1)
                     call simp_gen(intg,ff,pawrad(itypat))
                     intc(ll+1,iln,jln)=intg
                     have_intc(ll+1,iln,jln)=.true.
                   end if

                   do mm=-ll,ll
                     lm = LMC(ll,mm)
                     glmij=pawang%gntselect(lm,klm)
                     if (glmij>0) then
                       arg=ylmp(lm)*pawang%realgnt(glmij)*intc(ll+1,iln,jln)
                       radsumc(itypat)%value(1,ilmn,jlmn,i1,j2,i3,jkpt) = &
&                       radsumc(itypat)%value(1,ilmn,jlmn,i1,j2,i3,jkpt)+arg*real(ifac)
                       radsumc(itypat)%value(2,ilmn,jlmn,i1,j2,i3,jkpt) = &
&                       radsumc(itypat)%value(2,ilmn,jlmn,i1,j2,i3,jkpt)+arg*aimag(ifac)
                     end if
                   end do !mm
                 end do !ll
               end do !ilmn
             end do !jlmn

!            ===== Valence part =====
!            Need int1=\int phi_i phi_j jl (pr) dr
!            and  int2=\int tphi_i tphi_j jl (pr) dr

             have_rad(:,:)= .FALSE.;radint1=zero;radint2=zero;radint3=zero

             do klmn=1,pawtab(itypat)%lmn2_size
               klm=pawtab(itypat)%indklmn(1,klmn);kln=pawtab(itypat)%indklmn(2,klmn)
               llmin=pawtab(itypat)%indklmn(3,klmn);llmax=pawtab(itypat)%indklmn(4,klmn)

               do ll=llmin,llmax,2
                 ifac=ilfac(mod(ll,4))

                 if (.not.have_rad(ll+1,kln)) then
                   ff(1:mesh_size)=pawtab(itypat)%phiphj(1:mesh_size,kln)*j_bessel(1:mesh_size,ll+1)
                   call simp_gen(intg,ff,pawrad(itypat))
                   radint1(ll+1,kln)=intg
                   ff(1:mesh_size)=pawtab(itypat)%tphitphj(1:mesh_size,kln)*j_bessel(1:mesh_size,ll+1)
                   call simp_gen(intg,ff,pawrad(itypat))
                   radint2(ll+1,kln)=intg
                   ff(1:mesh_size)=(pawtab(itypat)%phiphj  (1:mesh_size,kln) &
&                   -pawtab(itypat)%tphitphj(1:mesh_size,kln))&
&                   *j_bessel(1:mesh_size,ll+1)
                   call simp_gen(intg,ff,pawrad(itypat))
                   radint3(ll+1,kln)=intg
                   have_rad(ll+1,kln)=.true.
                 end if

                 do mm=-ll,ll
                   lm = LMC(ll,mm)
                   glmij=pawang%gntselect(lm,klm)
                   if (glmij>0) then
                     arg=ylmp(lm)*pawang%realgnt(glmij)
                     radsum1(itypat)%value(1,klmn,i1,j2,i3,jkpt) = &
&                     radsum1(itypat)%value(1,klmn,i1,j2,i3,jkpt)+real(ifac) *arg*radint1(ll+1,kln)
                     radsum1(itypat)%value(2,klmn,i1,j2,i3,jkpt) = &
&                     radsum1(itypat)%value(2,klmn,i1,j2,i3,jkpt)+aimag(ifac)*arg*radint1(ll+1,kln)
                     radsum2(itypat)%value(1,klmn,i1,j2,i3,jkpt) = &
&                     radsum2(itypat)%value(1,klmn,i1,j2,i3,jkpt)+real(ifac) *arg*radint2(ll+1,kln)
                     radsum2(itypat)%value(2,klmn,i1,j2,i3,jkpt) = &
&                     radsum2(itypat)%value(2,klmn,i1,j2,i3,jkpt)+aimag(ifac)*arg*radint2(ll+1,kln)
                     radsum3(itypat)%value(1,klmn,i1,j2,i3,jkpt) = &
&                     radsum3(itypat)%value(1,klmn,i1,j2,i3,jkpt)+real(ifac) *arg*radint3(ll+1,kln)
                     radsum3(itypat)%value(2,klmn,i1,j2,i3,jkpt) = &
&                     radsum3(itypat)%value(2,klmn,i1,j2,i3,jkpt)+aimag(ifac)*arg*radint3(ll+1,kln)
                   end if
                 end do !mm
               end do !ll
             end do !klmn

           end do ! end loop over i1
         end do ! end loop over i2
       end do ! end loop over i3
     end do ! end loop over ikpt

     ABI_DEALLOCATE(ff)
     ABI_DEALLOCATE(jbes)

     ABI_DEALLOCATE(j_bessel)
     ABI_DEALLOCATE(ylmp)

     ABI_DEALLOCATE(intc)
     ABI_DEALLOCATE(have_intc)

     ABI_DEALLOCATE(radint1)
     ABI_DEALLOCATE(radint2)
     ABI_DEALLOCATE(radint3)
     ABI_DEALLOCATE(have_rad)

     call xmpi_sum(radsumc(itypat)%value,mpi_enreg%comm_band,ierr)
     call xmpi_sum(radsum1(itypat)%value,mpi_enreg%comm_band,ierr)
     call xmpi_sum(radsum2(itypat)%value,mpi_enreg%comm_band,ierr)
     call xmpi_sum(radsum3(itypat)%value,mpi_enreg%comm_band,ierr)

   end do ! end loop over atom types
 end if ! PAW

!Allocate main memory
 ABI_ALLOCATE(rho_contrib,(cplex*nfft))
 ABI_ALLOCATE(rho_contrib_g,(cplex,nfft))
 ABI_ALLOCATE(rho_contrib_paw1,(cplex,nfft))
 ABI_ALLOCATE(rho_contrib_paw2,(cplex,nfft))
 ABI_ALLOCATE(rho_contrib_paw3,(cplex,nfft))

 ABI_ALLOCATE(rho_moment_v1,(nfft,my_ngrid))
 ABI_ALLOCATE(rho_moment_v2,(nfft,my_ngrid))
 ABI_ALLOCATE(rho_moment_core,(nfft,my_ngrid))
 ABI_ALLOCATE(rho_pw,(nfft,my_ngrid))
 rho_moment_v1=zero;rho_moment_v2=zero
 rho_pw=zero;rho_moment_core=zero

!Prepare gamma(r) for the state independent scheme
 ABI_ALLOCATE(gamma,(nfft,2))
 if (.not.state_dependent) then
   ABI_ALLOCATE(rhoe,(nfft,1))
   ABI_ALLOCATE(rhop,(nfft,1))
   if (electronpositron%particle==EP_ELECTRON) then
     rhoe(:,1)=rhor_ep_(:,1);rhop(:,1)=rhor_(:,1)
   else if (electronpositron%particle==EP_POSITRON) then
     rhoe(:,1)=rhor_(:,1);rhop(:,1)=rhor_ep_(:,1)
   end if
   call mkdenpos(iwarn,nfft,1,1,rhoe(:,1),dtset%xc_denpos)
   call mkdenpos(iwarn,nfft,1,1,rhop(:,1),dtset%xc_denpos)
   call gammapositron_fft(electronpositron,gamma,Crystal%gprimd,igamma,mpi_enreg,&
&   n3xccc,nfft,ngfft,rhoe(:,1),rhop(:,1),xccc3d)
   ABI_DEALLOCATE(rhoe)
   ABI_DEALLOCATE(rhop)
 else
   gamma=one
 end if

!Some allocations for state-dependent scheme
 if (state_dependent) then
!  Fake MPI data to be used in poslifetime; allow only FFT parallelism
   call initmpi_seq(mpi_enreg_seq)
   mpi_enreg_seq%my_natom=dtset%natom
   call set_mpi_enreg_fft(mpi_enreg_seq,mpi_enreg%comm_fft,mpi_enreg%distribfft,&
&   mpi_enreg%me_g0,mpi_enreg%paral_kgb)
!  Allocate memory for state-dependent scheme
   ABI_ALLOCATE(rhor_dop_el,(nfft))
   if (dtset%usepaw==1) then
     ABI_DATATYPE_ALLOCATE(pawrhoij_dop_el,(dtset%natom))
     call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,nspden_rhoij=nspden_rhoij,&
&            nspden=dtset%nspden,spnorb=dtset%pawspnorb,cpxocc=dtset%pawcpxocc)
     call pawrhoij_alloc(pawrhoij_dop_el,cplex_rhoij,nspden_rhoij,&
                         dtset%nspinor,dtset%nsppol,dtset%typat,&
                         pawtab=pawtab,use_rhoij_=1,use_rhoijp=1)
!    Cancel distribution of PAW data over atomic sites
!    We use here pawrhoij because polifetime routine
!    detects by itself the particle described by pawrhoij
     if (mpi_enreg%my_natom<dtset%natom) then
       ABI_DATATYPE_ALLOCATE(pawrhoij_all,(dtset%natom))
       call pawrhoij_nullify(pawrhoij_all)
       call pawrhoij_gather(pawrhoij,pawrhoij_all,-1,mpi_enreg%comm_atom, &
&       with_rhoijres=.false.,with_rhoij_=.false.,with_lmnmix=.false.)
       ABI_DATATYPE_ALLOCATE(pawrhoij_ep_all,(dtset%natom))
       call pawrhoij_nullify(pawrhoij_ep_all)
       call pawrhoij_gather(electronpositron%pawrhoij_ep,pawrhoij_ep_all,-1,mpi_enreg%comm_atom, &
&       with_rhoijres=.false.,with_rhoij_=.false.,with_lmnmix=.false.)
     else
       pawrhoij_all => pawrhoij
       pawrhoij_ep_all => electronpositron%pawrhoij_ep
     end if
   end if
 end if

!==============================================================================
!================ Loop over positronic states =================================

!LOOP OVER k POINTS
 ibg_pos=0;icg_pos=0;ikg_pos=0;bdtot_index_pos=0;isppol_pos=1
 do ikpt_pos=1,merge(1,nkpt,kgamma_only_positron)

!  Extract data for this kpt_pos
   npw_k_pos=npwarr(ikpt_pos)
   wtk_k_pos=dtset%wtk(ikpt_pos); if (kgamma_only_positron) wtk_k_pos=one
   istwf_k_pos=dtset%istwfk(ikpt_pos)
   nband_k_pos=dtset%nband(ikpt_pos+(isppol_pos-1)*nkpt)
   nband_cprj_k_pos=nband_k_pos/nproc_band
   mykpt_pos=.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt_pos,1,nband_k_pos,&
&   isppol_pos,mpi_enreg%me_kpt))

!  Retrieve additional data for this kpt_pos
   ABI_ALLOCATE(occ_k_pos,(nband_k_pos))
   occ_k_pos(:)=occ_pos_ptr(1+bdtot_index_pos:nband_k_pos+bdtot_index_pos)
   nband_eff_pos=1
   do ib_pos=1,nband_k_pos
     if (occ_k_pos(ib_pos)>tol8) nband_eff_pos=ib_pos
   end do
   if (mod(nband_eff_pos,blocksize)/=0) nband_eff_pos=((nband_eff_pos/blocksize)+1)*blocksize

   nblock_band_eff_pos=nband_eff_pos/blocksize

   mcg_pos=npw_k_pos*my_nspinor*nband_eff_pos
   ABI_ALLOCATE(cg_k_pos,(2,mcg_pos))

   mcprj_k_pos=0
   if (dtset%usepaw==1) then
     nband_cprj_eff_pos=nband_eff_pos/nproc_band
     mcprj_k_pos=my_nspinor*nband_cprj_eff_pos
     ABI_DATATYPE_ALLOCATE(cprj_k_pos,(dtset%natom,mcprj_k_pos))
     call pawcprj_alloc(cprj_k_pos,0,dimcprj)
   end if

   if (mpi_enreg%paral_kgb==0) then
     ABI_ALLOCATE(gbound_pos,(2*dtset%mgfft+8,2))
     ABI_ALLOCATE(kg_k_pos,(3,npw_k_pos))
   else if (mykpt_pos) then
     nullify(bandfft_kpt_pos)
   else
     ABI_DATATYPE_ALLOCATE(bandfft_kpt_pos,)
     call bandfft_kpt_reset(bandfft_kpt_pos)
   end if

!  Exchange data (WF components) between procs
   if (mykpt_pos) then
     cg_k_pos(:,1:mcg_pos)=cg_pos_ptr(:,icg_pos+1:icg_pos+mcg_pos)
     if (mpi_enreg%paral_kgb==0) kg_k_pos(:,1:npw_k_pos)=kg(:,1+ikg_pos:npw_k_pos+ikg_pos)
     if (dtset%usepaw==1) then
       call pawcprj_get(Crystal%atindx1,cprj_k_pos,cprj_pos_ptr,dtset%natom,1,ibg_pos,ikpt_pos,iorder_cprj,&
&       isppol_pos,mband_cprj_pos,dtset%mkmem,dtset%natom,nband_cprj_eff_pos,nband_k_pos,my_nspinor,&
&       dtset%nsppol,dtfil%unpaw,mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
     end if
     if (mpi_enreg%paral_kgb/=0) then
       jj=mpi_enreg%my_kpttab(ikpt_pos)
       bandfft_kpt_pos => bandfft_kpt(jj)
     end if
     do ii=0,mpi_enreg%nproc_kpt-1
       if (ii/=mpi_enreg%me_kpt) then
         tag=ikpt_pos+(isppol_pos-1)*nkpt+2*nkpt*ii
         call xmpi_send(cg_k_pos,ii,tag,mpi_enreg%comm_kpt,ierr)
         tag=tag+nkpt*(1+2*mpi_enreg%nproc_kpt)
         if (mpi_enreg%paral_kgb==0) then
           call xmpi_send(kg_k_pos,ii,tag,mpi_enreg%comm_kpt,ierr)
         else
           call bandfft_kpt_mpi_send(bandfft_kpt_pos,ii,tag,mpi_enreg%comm_kpt,ierr,profile='fourwf')
         end if
         if (dtset%usepaw==1) then
           call pawcprj_mpi_send(dtset%natom,mcprj_k_pos,dimcprj,0,cprj_k_pos,ii,mpi_enreg%comm_kpt,ierr)
         end if
       end if
     end do
   else
     ii=0;if (allocated(mpi_enreg%proc_distrb)) ii=mpi_enreg%proc_distrb(ikpt_pos,1,isppol_pos)
     tag=ikpt_pos+(isppol_pos-1)*nkpt+2*nkpt*mpi_enreg%me_kpt
     call xmpi_recv(cg_k_pos,ii,tag,mpi_enreg%comm_kpt,ierr)
     tag=tag+nkpt*(1+2*mpi_enreg%nproc_kpt)
     if (mpi_enreg%paral_kgb==0) then
       call xmpi_recv(kg_k_pos,ii,tag,mpi_enreg%comm_kpt,ierr)
     else
       call bandfft_kpt_mpi_recv(bandfft_kpt_pos,ii,tag,mpi_enreg%comm_kpt,ierr)
     end if
     if (dtset%usepaw==1) then
       call pawcprj_mpi_recv(dtset%natom,mcprj_k_pos,dimcprj,0,cprj_k_pos,ii,mpi_enreg%comm_kpt,ierr)
     end if
   end if

   if (mpi_enreg%paral_kgb==0) then
     call sphereboundary(gbound_pos,istwf_k_pos,kg_k_pos,dtset%mgfft,npw_k_pos)
   end if

   ABI_ALLOCATE(cwaver_pos,(cplex*nfft))
   ABI_ALLOCATE(cwaver_pos_block,(cplex*nfft*bandpp))
   if (dtset%usepaw==1) then
     ABI_DATATYPE_ALLOCATE(cprj_pos,(dtset%natom,my_nspinor))
     call pawcprj_alloc(cprj_pos,0,dimcprj)
   end if

!  ============================================================================
!  Loops on positronic bands

   do iblock_pos=1,nblock_band_eff_pos
     ib_pos=1+(iblock_pos-1)*blocksize
     if (any(abs(occ_k_pos(ib_pos:ib_pos+blocksize-1))>tol8)) then

       ABI_ALLOCATE(cwaveg_pos,(2,npw_k_pos*blocksize))
       ABI_ALLOCATE(cwaveaug_pos,(2,n4,n5,n6*bandpp))
       ABI_ALLOCATE(denpot_dum,(n4,n5,n6))
       ABI_ALLOCATE(fofgout_dum,(2,npw_k_pos*blocksize))
       iwavef_pos=(iblock_pos-1)*npw_k_pos*blocksize
       cwaveg_pos(:,1:npw_k_pos*blocksize)= &
&       cg_k_pos(:,iwavef_pos+1:iwavef_pos+npw_k_pos*blocksize)

!      Get positronic wave function in real space
       option=0
       if (mpi_enreg%paral_kgb==0) then
         weight_pos=occ_k_pos(ib_pos)*wtk_k_pos
         call fourwf(1,denpot_dum,cwaveg_pos,fofgout_dum,cwaveaug_pos,&
&         gbound_pos,gbound_pos,istwf_k_pos,kg_k_pos,kg_k_pos,&
&         dtset%mgfft,mpi_enreg,1,ngfft,npw_k_pos,npw_k_pos,&
&         n4,n5,n6,option,tim_fourwf,weight_pos,weight_pos,&
&         use_gpu_cuda=dtset%use_gpu_cuda)
       else
         call prep_fourwf(denpot_dum,blocksize,cwaveg_pos,cwaveaug_pos,&
&         iblock_pos,istwf_k_pos,dtset%mgfft,mpi_enreg,nband_k_pos,&
&         bandpp,ngfft,npw_k_pos,n4,n5,n6,occ_k_pos,option,Crystal%ucvol,wtk_k_pos,&
&         bandfft_kpt_tab=bandfft_kpt_pos,use_gpu_cuda=dtset%use_gpu_cuda)
       end if

       cwaver_pos_block=zero
       do ii=1,bandpp
         j3=(ii-1)*n3
         indx0=1+(ii-1)*cplex*nfft
         do i3=1,n3
           if (me_fft==fftn3_distrib(i3)) then
             indx=indx0+cplex*n1*n2*(ffti3_local(i3)-1)
             do i2=1,n2
               do i1=1,n1
                 cwaver_pos_block(indx  )=cwaveaug_pos(1,i1,i2,i3+j3)
                 cwaver_pos_block(indx+1)=cwaveaug_pos(2,i1,i2,i3+j3)
                 indx=indx+2
               end do
             end do
           end if
         end do
       end do
       ABI_DEALLOCATE(fofgout_dum)
       ABI_DEALLOCATE(denpot_dum)
       ABI_DEALLOCATE(cwaveaug_pos)
       ABI_DEALLOCATE(cwaveg_pos)

!      At this stage, each band proc has bandpp bands in real space
!      (distributed on FFT procs)

!      ========================================================================
!      Compute core contribution for this positronic band (PAW only)

       if (dtset%usepaw==1) then
         do ibpp_pos=1,bandpp
           ib_cprj_pos=(iblock_pos-1)*bandpp+ibpp_pos
           weight_pos=occ_k_pos(ib_pos+ibpp_pos-1+me_band*bandpp)*wtk_k_pos
!       Calculate the annihilation rate for each core state for state dependent scheme
           iatm=0
           do itypat=1,dtset%ntypat
             mesh_size = pawtab(itypat)%mesh_size
             do iat=1,Crystal%nattyp(itypat)
               iatm=iatm+1;iatom=Crystal%atindx1(iatm)
               ABI_ALLOCATE(gammastate_c(iatom)%value,(lmncmax(itypat)))
               do jlmn=1,lmncmax(itypat)
                 jln = indlmncor(itypat)%value(5,jlmn)
                 contrib(:)=zero
                 ABI_ALLOCATE(rhocorej,(mesh_size))
                 rhocorej(1:mesh_size)=2*phicor(itypat)%value(1:mesh_size,jln)**2
                 call posratecore(dtset,electronpositron,iatom,dtset%natom,mesh_size,mpi_enreg_seq,&
&                 1,pawang,pawrad,pawrhoij_all,pawrhoij_ep_all,pawtab,ratec,rhocorej)

                 call posratecore(dtset,electronpositron,iatom,dtset%natom,mesh_size,mpi_enreg_seq,&
&                 2,pawang,pawrad,pawrhoij_all,pawrhoij_ep_all,pawtab,ratec_ipm,rhocorej)

                 gammastate_c(iatom)%value(jlmn)=ratec/ratec_ipm
                 ABI_DEALLOCATE(rhocorej)
               end do
             end do
           end do
           jkpt=0
           do ikpt=1,nkpt
             if (my_gridtab(ikpt)==0) cycle
             jkpt=jkpt+1
             do i3=1,n3
               ig3=i3-(i3/id3)*n3-1
               do i2=1,n2
                 if (me_fft==fftn2_distrib(i2)) then
                   j2=ffti2_local(i2)
                   ig2=i2-(i2/id2)*n2-1
                   indx=n1*(my_n2*(i3-1)+(j2-1))
                   do i1=1,n1
                     ig1=i1-(i1/id1)*n1-1
                     indx=indx+1

!                    Loop on atoms (type sorted)
                     iatm=0
                     do itypat=1,dtset%ntypat
                       lmn_size = pawtab(itypat)%lmn_size

                       do iat=1,Crystal%nattyp(itypat)
                         iatm=iatm+1;iatom=Crystal%atindx1(iatm)

                         pcart(:)=Crystal%gprimd(:,1)*real(ig1+dtset%kpt(1,ikpt))+&
&                         Crystal%gprimd(:,2)*real(ig2+dtset%kpt(2,ikpt))+&
&                         Crystal%gprimd(:,3)*real(ig3+dtset%kpt(3,ikpt))
                         pnorm=dsqrt(dot_product(pcart,pcart))
                         pr=dot_product(pcart,Crystal%xcart(:,iatom))
                         expipr(1)= cos(two_pi*pr)
                         expipr(2)=-sin(two_pi*pr)

!                        Loop on ij states
                         do jlmn = 1,lmncmax(itypat)
                           contrib(:)=zero
                           do ilmn = 1,lmn_size
                             radsumnfftc(1)=expipr(1)*radsumc(itypat)%value(1,ilmn,jlmn,i1,j2,i3,jkpt)&
&                             -expipr(2)*radsumc(itypat)%value(2,ilmn,jlmn,i1,j2,i3,jkpt)
                             radsumnfftc(2)=expipr(1)*radsumc(itypat)%value(2,ilmn,jlmn,i1,j2,i3,jkpt)&
&                             +expipr(2)*radsumc(itypat)%value(1,ilmn,jlmn,i1,j2,i3,jkpt)
                             cp_pos(:)=cprj_k_pos(iatom,ib_cprj_pos)%cp(:,ilmn)
                             contrib(1)=contrib(1)+four_pi*(cp_pos(1)*radsumnfftc(1) &
&                             -cp_pos(2)*radsumnfftc(2))
                             contrib(2)=contrib(2)+four_pi*(cp_pos(1)*radsumnfftc(2) &
&                             +cp_pos(2)*radsumnfftc(1))
                           end do ! end loop over ilmn
                           ! 2 - electron state weight for 2 spins
                           rho_moment_core(indx,jkpt) = rho_moment_core(indx,jkpt) &
&                           +gammastate_c(iatom)%value(jlmn)*2*weight_pos*(contrib(1)**2+contrib(2)**2)
                         end do ! end loop over jlmn

                       end do !end loop over atoms
                     end do !end loop over atom types

                   end do ! end loop over i1
                 end if ! end loop over i2
               end do
             end do ! end loop over i3
           end do ! jkpt
         end do ! ibpp_pos
       end if

!      We now loop over positronic bands inside a block
!      and select occupied ones
       do ibpp_pos=1,blocksize
         ib_pos=(iblock_pos-1)*blocksize+ibpp_pos
         occ_pos=occ_k_pos(ib_pos)
         if (abs(occ_pos)>tol8) then

!          Parallelism: dirty trick (broadcast bands) but there should be few positronic bands (~1)
           if (nproc_band>1) then
             iproc=(ibpp_pos-1)/bandpp
             if (me_band==iproc) then
               indx=mod((ibpp_pos-1),bandpp)*cplex*nfft
               cwaver_pos(1:cplex*nfft)=cwaver_pos_block(indx+1:indx+cplex*nfft)
             end if
             call xmpi_bcast(cwaver_pos,iproc,mpi_enreg%comm_band,ierr)
             if (dtset%usepaw==1) then
               if (me_band==iproc) then
                 indx=mod((ibpp_pos-1),bandpp)*my_nspinor
                 call pawcprj_copy(cprj_k_pos(:,indx+1:indx+my_nspinor),cprj_pos)
               end if
               call pawcprj_bcast(cprj_pos,dtset%natom,my_nspinor,dimcprj,0,iproc,&
&               mpi_enreg%comm_band,ierr)
             end if
           else
             cwaver_pos(1:cplex*nfft)=cwaver_pos_block(1:cplex*nfft)
             if (dtset%usepaw==1) then
               call pawcprj_copy(cprj_k_pos(:,(ib_pos-1)*my_nspinor+1:ib_pos*my_nspinor),cprj_pos)
             end if
           end if

!      ========================================================================
!      ================ Loop over electronic states ===========================

!          Loop over spins
           ibg=0;icg=0;ikg=0;bdtot_index=0
           do isppol=1,dtset%nsppol
!            Loop over k points
             ikg=0;jkpt=0
             do ikpt=1,nkpt

!              Extract data for this kpt_pos
               npw_k=npwarr(ikpt)
               wtk_k=dtset%wtk(ikpt)
               istwf_k=dtset%istwfk(ikpt)
               nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
               nband_cprj_k=nband_k/nproc_band
               mykpt=.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,&
&               isppol,mpi_enreg%me_kpt))

!              Select k-points for current proc
               if (mykpt) then

!                Retrieve additional data for this kpt_pos
                 jkpt=jkpt+1
                 ABI_ALLOCATE(occ_k,(nband_k))
                 occ_k(:)=occ_ptr(1+bdtot_index:nband_k+bdtot_index)

                 mcprj_k=0
                 if (dtset%usepaw==1) then
                   mcprj_k=my_nspinor*nband_cprj_k
                   ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,mcprj_k))
                   call pawcprj_alloc(cprj_k,0,dimcprj)
                   call pawcprj_get(Crystal%atindx1,cprj_k,cprj_ptr,dtset%natom,1,ibg,ikpt,iorder_cprj,&
&                   isppol,mband_cprj,dtset%mkmem,dtset%natom,nband_cprj_k,nband_cprj_k,my_nspinor,&
&                   dtset%nsppol,dtfil%unpaw,mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
                 end if

                 if (mpi_enreg%paral_kgb==0) then
                   ABI_ALLOCATE(gbound,(2*dtset%mgfft+8,2))
                   ABI_ALLOCATE(kg_k,(3,npw_k))
                   kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
                   call sphereboundary(gbound,istwf_k,kg_k,dtset%mgfft,npw_k)
                 else
                   jj=mpi_enreg%my_kpttab(ikpt)
                   bandfft_kpt_el => bandfft_kpt(jj)
                 end if

                 ABI_ALLOCATE(cwaver,(cplex*nfft*bandpp))

!                ==================================================================
!                Loops on electronic bands

                 do iblock=1,nblock_band
                   ib=1+(iblock-1)*blocksize

                   if (any(abs(occ_k(ib:ib+blocksize-1))>tol8)) then

!                    Retrieve electronic wave function
                     ABI_ALLOCATE(cwaveg,(2,npw_k*blocksize))
                     ABI_ALLOCATE(cwaveaug,(2,n4,n5,n6*bandpp))
                     ABI_ALLOCATE(denpot_dum,(n4,n5,n6))
                     ABI_ALLOCATE(fofgout_dum,(2,npw_k*blocksize))
                     iwavef=(iblock-1)*npw_k*blocksize
                     cwaveg(:,1:npw_k*blocksize)= &
&                     cg_ptr(:,icg+iwavef+1:icg+iwavef+npw_k*blocksize)

!                    Get electronic wave function in real space
                     option=0
                     if (mpi_enreg%paral_kgb==0) then
                       weight=occ_k(ib)*wtk_k
                       call fourwf(1,denpot_dum,cwaveg,fofgout_dum,cwaveaug,&
&                       gbound,gbound,istwf_k,kg_k,kg_k,&
&                       dtset%mgfft,mpi_enreg,1,ngfft,npw_k,npw_k,&
&                       n4,n5,n6,option,tim_fourwf,weight,weight,&
&                       use_gpu_cuda=dtset%use_gpu_cuda)
                     else
                       call prep_fourwf(denpot_dum,blocksize,cwaveg,cwaveaug,&
&                       iblock,istwf_k,dtset%mgfft,mpi_enreg,nband_k,&
&                       bandpp,ngfft,npw_k,n4,n5,n6,occ_k,option,Crystal%ucvol,wtk_k,&
&                       bandfft_kpt_tab=bandfft_kpt_el,use_gpu_cuda=dtset%use_gpu_cuda)
                     end if

                     cwaver=zero
                     do ii=1,bandpp
                       j3=(ii-1)*n3
                       indx0=1+(ii-1)*cplex*nfft
                       do i3=1,n3
                         if (me_fft==fftn3_distrib(i3)) then
                           indx=indx0+cplex*n1*n2*(ffti3_local(i3)-1)
                           do i2=1,n2
                             do i1=1,n1
                               cwaver(indx  )=cwaveaug(1,i1,i2,i3+j3)
                               cwaver(indx+1)=cwaveaug(2,i1,i2,i3+j3)
                               indx=indx+2
                             end do
                           end do
                         end if
                       end do
                     end do
                     ABI_DEALLOCATE(fofgout_dum)
                     ABI_DEALLOCATE(denpot_dum)
                     ABI_DEALLOCATE(cwaveaug)
                     ABI_DEALLOCATE(cwaveg)
!                    At this stage, each band proc has bandpp bands in real space
!                   (distributed on FFT procs)

!                    We now loop on the bandpp bands
!                    and select occupied ones
                     do ibpp=1,bandpp
                       occ_el=occ_k(ib+ibpp-1+me_band*bandpp)
                       if (abs(occ_el)>tol8) then

!                        ==============================================================
!                        Compute state-dependent annihilation rate
!                        Avoid parallelism over kpt/bands/atoms
                         gammastate=one;rate_paw=one
                         if (state_dependent) then
                           weight=occ_el*wtk_k
                           ib_cprj=(iblock-1)*bandpp+ibpp
                           indx=1+(ibpp-1)*cplex*nfft
                           do ii=1,nfft
                             rhor_dop_el(ii)=weight*(cwaver(indx)*cwaver(indx)+cwaver(indx+1)*cwaver(indx+1))
                             indx=indx+2
                           end do
                           if (dtset%usepaw==1) then
                             do iatom=1,dtset%natom
                               pawrhoij_dop_el(iatom)%rhoij_=zero
                             end do
                             cplex_rhoij=2;if (istwf_k>1) cplex_rhoij=1
                             usetimerev=(dtset%kptopt>0.and.dtset%kptopt<3)
                             call pawaccrhoij(Crystal%atindx,cplex_rhoij,cprj_k(:,ib_cprj),cprj_k(:,ib_cprj),0,isppol,&
&                             dtset%natom,dtset%natom,dtset%nspinor,occ_el,1,pawrhoij_dop_el,usetimerev,wtk_k)
!                            Is it correct to apply symetries here (on a single band)?
!                            If not, call pawrhoij_symrhoij with nsym=1
                             call pawrhoij_symrhoij(pawrhoij_dop_el,pawrhoij_dop_el,1,Crystal%gprimd,&
&                             Crystal%indsym,0,dtset%natom,Crystal%nsym,dtset%ntypat,1,pawang,-10001,&
&                             pawtab,Crystal%rprimd,Crystal%symafm,Crystal%symrec,dtset%typat)
                           end if
!                          Has to call poslifetime in sequential because we are in a parallel section
!                          Only FFT parallelism is allowed
                           call poslifetime(dtset,electronpositron,Crystal%gprimd,dtset%natom,mpi_enreg_seq,n3xccc,&
&                           nfft,ngfft,nhat,2,pawang,pawrad,pawrhoij_all,pawtab,rate,rate_paw,rhor,Crystal%ucvol,xccc3d,&
&                           rhor_dop_el=rhor_dop_el,pawrhoij_dop_el=pawrhoij_dop_el,pawrhoij_ep=pawrhoij_ep_all)
                           call poslifetime(dtset,electronpositron,Crystal%gprimd,dtset%natom,mpi_enreg_seq,n3xccc,&
&                           nfft,ngfft,nhat,3,pawang,pawrad,pawrhoij_all,pawtab,rate_ipm,rate_paw_ipm,rhor,Crystal%ucvol,xccc3d,&
&                           rhor_dop_el=rhor_dop_el,pawrhoij_dop_el=pawrhoij_dop_el,pawrhoij_ep=pawrhoij_ep_all)
                           gammastate=rate/rate_ipm
                           rate_paw=rate_paw/rate_paw_ipm
                         end if

!                        ==============================================================
!                        Compute plane-wave contribution to momentum distribution

!                        Compute Psi^+(r) * Psi^-(r) * gamma(r) in real space
                         rho_contrib(:)=zero
                         indx=(ibpp-1)*cplex*nfft
                         if (cplex==2) then
                           do jj=1,nfft
                             ii=2*jj-1
                             rho_contrib(ii)  =sqrt(gamma(jj,2))*(cwaver_pos(ii)*cwaver(indx+ii)&
&                             -wf_fact*cwaver_pos(ii+1)*cwaver(indx+ii+1))
                             rho_contrib(ii+1)=sqrt(gamma(jj,2))*(cwaver_pos(ii)*cwaver(indx+ii+1) &
&                             +wf_fact*cwaver_pos(ii+1)*cwaver(indx+ii))
                           end do
                         else
                           do ii=1,nfft
                             rho_contrib(ii)=sqrt(gamma(ii,2))*cwaver_pos(ii)*cwaver(indx+ii)
                           end do
                         end if

!                        FFT of (Psi+.Psi-.gamma) to get Intg[(Psi+.Psi-.gamma).exp(-igr)]
                         call fourdp(cplex,rho_contrib_g,rho_contrib,-1,mpi_enreg,nfft,1,ngfft,&
&                         tim_fourdp)

                         rho_pw(1:nfft,jkpt)=rho_pw(1:nfft,jkpt) +gammastate*occ_el*occ_pos &
&                         *(rho_contrib_g(1,1:nfft)**2+rho_contrib_g(2,1:nfft)**2)

!                        ==============================================================
!                        Compute PAW on-site contribution to momentum distribution

                         if (dtset%usepaw==1) then

                           rho_contrib_paw1(:,:)= zero
                           rho_contrib_paw2(:,:)= zero
                           rho_contrib_paw3(:,:)= zero

                           ib_cprj=(iblock-1)*bandpp+ibpp

!                          Loop on moments
                           indx=0
                           do i3=1,n3
                             ig3=i3-(i3/id3)*n3-1
                             do i2=1,n2
                               if (me_fft==fftn2_distrib(i2)) then
                                 j2=ffti2_local(i2)
                                 ig2=i2-(i2/id2)*n2-1
                                 indx=n1*(my_n2*(i3-1)+(j2-1))
                                 do i1=1,n1
                                   ig1=i1-(i1/id1)*n1-1
                                   indx=indx+1

                                   pcart(:)=Crystal%gprimd(:,1)*real(ig1+dtset%kpt(1,ikpt))+&
&                                   Crystal%gprimd(:,2)*real(ig2+dtset%kpt(2,ikpt))+&
&                                   Crystal%gprimd(:,3)*real(ig3+dtset%kpt(3,ikpt))
                                   pnorm=dsqrt(dot_product(pcart,pcart))

!                                  Loop on atoms (type-sorted)
                                   iatm=0
                                   do itypat=1,dtset%ntypat
                                     lmn_size=pawtab(itypat)%lmn_size
                                     lmn2_size=pawtab(itypat)%lmn2_size
                                     ABI_ALLOCATE(radsumnfft1,(2,lmn2_size))
                                     ABI_ALLOCATE(radsumnfft2,(2,lmn2_size))
                                     ABI_ALLOCATE(radsumnfft3,(2,lmn2_size))

                                     do iat=1,Crystal%nattyp(itypat)
                                       iatm=iatm+1;iatom=Crystal%atindx1(iatm)

                                       pr=dot_product(pcart,Crystal%xcart(:,iatom))
                                       expipr(1)= cos(two_pi*pr)
                                       expipr(2)=-sin(two_pi*pr)

                                       do klmn=1,lmn2_size
                                         radsumnfft1(1,klmn)=expipr(1)*radsum1(itypat)%value(1,klmn,i1,j2,i3,jkpt)&
&                                         -expipr(2)*radsum1(itypat)%value(2,klmn,i1,j2,i3,jkpt)
                                         radsumnfft1(2,klmn)=expipr(1)*radsum1(itypat)%value(2,klmn,i1,j2,i3,jkpt)&
&                                         +expipr(2)*radsum1(itypat)%value(1,klmn,i1,j2,i3,jkpt)
                                         radsumnfft2(1,klmn)=expipr(1)*radsum2(itypat)%value(1,klmn,i1,j2,i3,jkpt)&
&                                         -expipr(2)*radsum2(itypat)%value(2,klmn,i1,j2,i3,jkpt)
                                         radsumnfft2(2,klmn)=expipr(1)*radsum2(itypat)%value(2,klmn,i1,j2,i3,jkpt)&
&                                         +expipr(2)*radsum2(itypat)%value(1,klmn,i1,j2,i3,jkpt)
                                         radsumnfft3(1,klmn)=expipr(1)*radsum3(itypat)%value(1,klmn,i1,j2,i3,jkpt)&
&                                         -expipr(2)*radsum3(itypat)%value(2,klmn,i1,j2,i3,jkpt)
                                         radsumnfft3(2,klmn)=expipr(1)*radsum3(itypat)%value(2,klmn,i1,j2,i3,jkpt)&
&                                         +expipr(2)*radsum3(itypat)%value(1,klmn,i1,j2,i3,jkpt)
                                       end do

!                                      Loop on ij states
                                       do ilmn = 1, lmn_size
                                         i0lmn = ilmn*(ilmn-1)/2
                                         do jlmn = 1, lmn_size
                                           klmn = i0lmn+jlmn
                                           if (jlmn>ilmn) then
                                             i0lmn=jlmn*(jlmn-1)/2; klmn=i0lmn+ilmn
                                           end if
!                                          Transform 3-dimentional radsum to 1-dimentional radsumnfft
                                           cp(:)=cprj_k(iatom,ib_cprj)%cp(:,ilmn)
                                           cp_pos(:)=cprj_pos(iatom,1)%cp(:,jlmn)
                                           cp11= cp(1)*cp_pos(1)
                                           cp22= cp(2)*cp_pos(2)*wf_fact
                                           cp12= cp(1)*cp_pos(2)*wf_fact
                                           cp21= cp(2)*cp_pos(1)
                                           cpr=cp11-cp22 ; cpi=cp12+cp21
                                           rho_contrib_paw1(1,indx) = rho_contrib_paw1(1,indx) &
&                                           + four_pi*(cpr*radsumnfft1(1,klmn)-cpi*radsumnfft1(2,klmn))
                                           rho_contrib_paw1(2,indx) = rho_contrib_paw1(2,indx) &
&                                           + four_pi*(cpr*radsumnfft1(2,klmn)+cpi*radsumnfft1(1,klmn))
                                           rho_contrib_paw2(1,indx) = rho_contrib_paw2(1,indx) &
&                                           + four_pi*(cpr*radsumnfft2(1,klmn)-cpi*radsumnfft2(2,klmn))
                                           rho_contrib_paw2(2,indx) = rho_contrib_paw2(2,indx) &
&                                           + four_pi*(cpr*radsumnfft2(2,klmn)+cpi*radsumnfft2(1,klmn))
                                           rho_contrib_paw3(1,indx) = rho_contrib_paw3(1,indx) &
&                                           + four_pi*(cpr*radsumnfft3(1,klmn)-cpi*radsumnfft3(2,klmn))
                                           rho_contrib_paw3(2,indx) = rho_contrib_paw3(2,indx) &
&                                           + four_pi*(cpr*radsumnfft3(2,klmn)+cpi*radsumnfft3(1,klmn))
                                         end do ! end loop over jlmn
                                       end do ! end loop over ilmn

                                     end do !end loop over atoms

                                     ABI_DEALLOCATE(radsumnfft1)
                                     ABI_DEALLOCATE(radsumnfft2)
                                     ABI_DEALLOCATE(radsumnfft3)
                                   end do !end loop over atom types

                                   rho_moment_v1(indx,jkpt) = rho_moment_v1(indx,jkpt) &
&                                   +occ_el*occ_pos &
&                                   *(gammastate*(rho_contrib_g(1,indx)**2+rho_contrib_g(2,indx)**2) &
&                                   +rate_paw*(rho_contrib_paw1(1,indx)**2+rho_contrib_paw1(2,indx)**2 &
&                                   -rho_contrib_paw2(1,indx)**2-rho_contrib_paw2(2,indx)**2))
                                   rho_moment_v2(indx,jkpt) = rho_moment_v2(indx,jkpt) &
&                                   +occ_el*occ_pos*gammastate &
&                                   *((rho_contrib_g(1,indx)+rho_contrib_paw3(1,indx))**2+&
&                                   (rho_contrib_g(2,indx)+rho_contrib_paw3(2,indx))**2)

                                 end do ! end loop over i1

                               end if ! end loop over i2
                             end do
                           end do ! end loop over i3

                         end if ! PAW

!                      ================================================================
!                      End loops on electronic bands

                       end if ! occ>1.e-8
                     end do ! ibpp
                   end if ! occ_block>1.e-8
                 end do ! iblock

!                End loops over k points and spins (electrons)
                 icg = icg + npw_k*my_nspinor*nband_k
                 ibg = ibg + my_nspinor*nband_cprj_k
                 ikg = ikg + npw_k

                 ABI_DEALLOCATE(cwaver)
                 ABI_DEALLOCATE(occ_k)
                 if (mpi_enreg%paral_kgb==0) then
                   ABI_DEALLOCATE(kg_k)
                   ABI_DEALLOCATE(gbound)
                 else
                   nullify(bandfft_kpt_el)
                 end if
                 if (dtset%usepaw==1) then
                   call pawcprj_free(cprj_k)
                   ABI_DATATYPE_DEALLOCATE(cprj_k)
                 end if

               end if ! mykpt
               bdtot_index=bdtot_index+nband_k
             end do ! ikpt
           end do ! isppol

!          ================================================================
!          End loops on positronic bands

         end if ! occ>1.e-8
       end do ! ibpp_pos
     end if ! occ(block)>1.e-8
   end do ! iblock_pos

!  End loop over k points (positron)
   if (mykpt_pos) then
     icg_pos = icg_pos + npw_k_pos*my_nspinor*nband_k_pos
     ibg_pos = ibg_pos + my_nspinor*nband_cprj_k_pos
     ikg_pos = ikg_pos + npw_k_pos
   end if
   bdtot_index_pos=bdtot_index_pos+nband_k_pos

   ABI_DEALLOCATE(cwaver_pos)
   ABI_DEALLOCATE(cwaver_pos_block)
   ABI_DEALLOCATE(cg_k_pos)
   ABI_DEALLOCATE(occ_k_pos)
   if (mpi_enreg%paral_kgb==0) then
     ABI_DEALLOCATE(kg_k_pos)
     ABI_DEALLOCATE(gbound_pos)
   else if (mykpt_pos) then
     nullify(bandfft_kpt_pos)
   else
     call bandfft_kpt_destroy(bandfft_kpt_pos)
     ABI_DATATYPE_DEALLOCATE(bandfft_kpt_pos)
   end if
   if (dtset%usepaw==1) then
     call pawcprj_free(cprj_pos)
     ABI_DATATYPE_DEALLOCATE(cprj_pos)
     call pawcprj_free(cprj_k_pos)
     ABI_DATATYPE_DEALLOCATE(cprj_k_pos)
   end if

 end do ! ikpt_pos

!================================================================
!Final computations and printing

!In case of parallelism, sum over the communicator(s)
 if (nproc_band>1) then
   ABI_ALLOCATE(mpibuf,(3*nfft,my_ngrid))
   do jkpt=1,my_ngrid
     mpibuf(       1:  nfft,jkpt)=rho_moment_v1(1:nfft,jkpt)
     mpibuf(  nfft+1:2*nfft,jkpt)=rho_moment_v2(1:nfft,jkpt)
     mpibuf(2*nfft+1:3*nfft,jkpt)=rho_pw       (1:nfft,jkpt)
   end do
   call xmpi_sum(mpibuf,mpi_enreg%comm_band,ierr)
   do jkpt=1,my_ngrid
     rho_moment_v1(1:nfft,jkpt)=mpibuf(       1:  nfft,jkpt)
     rho_moment_v2(1:nfft,jkpt)=mpibuf(  nfft+1:2*nfft,jkpt)
     rho_pw(1:nfft,jkpt)       =mpibuf(2*nfft+1:3*nfft,jkpt)
   end do
   ABI_DEALLOCATE(mpibuf)
 end if
 if (dtset%usepaw==1) then
   call xmpi_sum(rho_moment_core,mpi_enreg%comm_band,ierr)
 end if

!Add valence and core contributions
 if (dtset%usepaw==1) then
   if (dtset%nsppol==2.and.my_nsppol==1) rho_moment_core(:,:)=half*rho_moment_core(:,:)
   rho_moment_v1(:,:)=rho_moment_v1(:,:)+rho_moment_core(:,:)
   rho_moment_v2(:,:)=rho_moment_v2(:,:)+rho_moment_core(:,:)
 end if

 units_=pi*(one/InvFineStruct)**3/Time_Sec/1.e12_dp/electronpositron%posocc
 scale_=(two_pi**2)/(Crystal%ucvol**two_thirds)

!Integrate rho_moment over p
 buf(1)=sum(rho_moment_v1(1:nfft,1:my_ngrid))
 buf(2)=sum(rho_moment_v2(1:nfft,1:my_ngrid))
 buf(3)=sum(rho_moment_core(1:nfft,1:my_ngrid))
 buf(4)=sum(rho_pw(1:nfft,1:my_ngrid))
 call xmpi_sum(buf,mpi_enreg%comm_kpt,ierr)
 call xmpi_sum(buf,mpi_enreg%comm_fft,ierr)
 lambda_v1=buf(1)*units_/Crystal%ucvol/nkpt
 lambda_v2=buf(2)*units_/Crystal%ucvol/nkpt
 lambda_core=buf(3)*units_/Crystal%ucvol/nkpt
 lambda_pw=buf(4)*units_/Crystal%ucvol/nkpt

!Write result in _DOPPLER file
!Requires MPI-IO if nproc_fft>1
 if (me_band==0) then
   if (me_kpt==0) then
     filename_dop=trim(dtfil%filnam_ds(4))//'_DOPPLER'
     vec=sqrt(dot_product(Crystal%gprimd(:,3),Crystal%gprimd(:,3)))
     ABI_ALLOCATE(pcart_k,(3,nfft))
     ABI_ALLOCATE(rho_moment_k,(nfft))
     if (dtset%nsppol==2) then
       ABI_ALLOCATE(rho_moment_k2,(nfft))
     end if
     if (accessfil==IO_MODE_FORTRAN) then  ! >>>>> Fortran access
!      Open file and write first line
       ierr=open_file(filename_dop,msg,newunit=unit_doppler,form='unformatted')
       write(unit_doppler) nfft,nkpt,Crystal%ucvol,Crystal%rprimd(:,:)
     else                                 ! >>>>> MPI-IO access
       unit_doppler=get_unit()
!      Open file and write first line
       call WffOpen(IO_MODE_MPI,mpi_enreg%comm_fft,filename_dop,ierr,wff,0,me_fft,unit_doppler)
       if (me_fft==0) then
         call xderiveWRecInit(wff,ierr)
         call xderiveWrite(wff,n1*n2*n3,ierr)
         call xderiveWrite(wff,nkpt,ierr)
         call xderiveWrite(wff,Crystal%ucvol,ierr)
         call xderiveWrite(wff,Crystal%rprimd(:,:),ierr)
         call xderiveWRecEnd(wff,ierr)
       else
         call xmoveOff(wff,n_int=2,n_dp=10,n_mark=2)
       end if
!      Store table of FFT points treated by current proc
       ABI_ALLOCATE(my_ffttab,(nfft))
       my_ffttab=0
       do i3=1,n3
         do i2=1,n2
           if (me_fft==fftn2_distrib(i2)) then
             indx0=n1*(n2*(i3-1)+(i2-1))
             indx=n1*(my_n2*(i3-1)+(ffti2_local(i2)-1))
             my_ffttab(indx+1:indx+n1)=(/(indx0+ii,ii=1,n1)/)
           end if
         end do
       end do
       ABI_ALLOCATE(mpibuf,(1,nfft))
     end if
   end if

   jkpt=0
   do ikpt=1,nkpt
     if (nproc_kpt==1) then
       rho_moment_k(1:nfft)=rho_moment_v2(1:nfft,ikpt)
     else
       if (my_gridtab(ikpt)/=0) jkpt=jkpt+1
       if (me_kpt==0) then
         if (my_gridtab(ikpt)==0) then
           tag=ikpt;iproc=mpi_enreg%proc_distrb(ikpt,1,1)
           call xmpi_recv(rho_moment_k,iproc,tag,mpi_enreg%comm_kpt,ierr)
           if (dtset%nsppol==2) then
             tag=2*ikpt;iproc=mpi_enreg%proc_distrb(ikpt,1,2)
             call xmpi_recv(rho_moment_k2,iproc,tag,mpi_enreg%comm_kpt,ierr)
             rho_moment_k(1:nfft)=rho_moment_k(1:nfft)+rho_moment_k2(1:nfft)
           end if
         else if (any(mpi_enreg%my_isppoltab(:)==1)) then
           rho_moment_k(1:nfft)=rho_moment_v2(1:nfft,jkpt)
           if (dtset%nsppol==2) then
             ii=2;if (mpi_enreg%my_isppoltab(2)==1) ii=1
             tag=ii*ikpt;iproc=mpi_enreg%proc_distrb(ikpt,1,ii)
             call xmpi_recv(rho_moment_k2,iproc,tag,mpi_enreg%comm_kpt,ierr)
             rho_moment_k(1:nfft)=rho_moment_k(1:nfft)+rho_moment_k2(1:nfft)
           end if
         end if
       else if (my_gridtab(ikpt)/=0) then
         if (mpi_enreg%my_isppoltab(1)==1) then
           tag=ikpt
           call xmpi_send(rho_moment_v2(1:nfft,jkpt),0,tag,mpi_enreg%comm_kpt,ierr)
         end if
         if (dtset%nsppol==2.and.mpi_enreg%my_isppoltab(2)==1) then
           tag=2*ikpt
           call xmpi_send(rho_moment_v2(1:nfft,jkpt),0,tag,mpi_enreg%comm_kpt,ierr)
         end if
       end if
     end if ! nproc_kpt>1
     if (me_kpt==0) then
       indx=0
       do i3=1,n3
         ig3=i3-(i3/id3)*n3-1
         do i2=1,n2
           if (me_fft/=fftn2_distrib(i2)) cycle
           ig2=i2-(i2/id2)*n2-1
           do i1=1,n1
             ig1=i1-(i1/id1)*n1-1
             indx=indx+1
             pcart_k(:,indx)=Crystal%gprimd(:,1)*real(ig1+dtset%kpt(1,ikpt))+&
&             Crystal%gprimd(:,2)*real(ig2+dtset%kpt(2,ikpt))+&
&             Crystal%gprimd(:,3)*real(ig3+dtset%kpt(3,ikpt))
           end do
         end do
       end do
       if (accessfil==IO_MODE_FORTRAN) then
         write(unit_doppler) pcart_k(1:3,1:nfft),rho_moment_k(1:nfft)
       else
         mpibuf(1,1:nfft)=rho_moment_k(1:nfft)
         call xderiveWRecInit(wff,ierr)
         call xderiveWrite(wff,pcart_k,3,nfft,mpi_enreg%comm_fft,my_ffttab,ierr)
         call xderiveWrite(wff,mpibuf ,1,nfft,mpi_enreg%comm_fft,my_ffttab,ierr)
         call xderiveWRecEnd(wff,ierr)
       end if
     end if
   end do
   if (me_kpt==0) then
     ABI_DEALLOCATE(pcart_k)
     ABI_DEALLOCATE(rho_moment_k)
     if (dtset%nsppol==2) then
       ABI_DEALLOCATE(rho_moment_k2)
     end if
     if (accessfil==IO_MODE_FORTRAN) then
       ierr=close_unit(unit_doppler,msg)
     else
       call WffClose(wff,ierr)
       ABI_DEALLOCATE(my_ffttab)
       ABI_DEALLOCATE(mpibuf)
     end if
   end if
 end if ! me_band==0

!Write results
 write(msg,'(7a)') &
& ' Computation of electron-positron pairs momentum distribution completed.',ch10,&
& '-File ',trim(filename_dop),' has been created.',ch10,&
& '-Use ~abinit/scripts/post_processing/posdopspectra.F90 to process it.'
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')
 msg=' Some e-p annihilation rates (ns-1) obtained by integration of e-p pairs momentum distribution:'
 call wrtout(std_out,msg,'COLL')
 write(msg,'(a,es22.12,3(2a,es22.12))') &
& '   Lambda (from module of sum of PAW contrib.)  = ',lambda_v2*1000._dp,ch10,&
& '     = lambda_core: ',lambda_core*1000._dp,ch10,&
& '      +lambda_pw  : ',lambda_pw*1000._dp,ch10,&
& '      +lambda_paw : ',(lambda_v2-lambda_core-lambda_pw)*1000._dp
 call wrtout(std_out,msg,'COLL')
 write(msg,'(4(a,es22.12,a))') &
& '   Lambda (from sum of modules of PAW contrib.) = ',lambda_v1*1000._dp,ch10,&
& '     = lambda_core: ',lambda_core*1000._dp,ch10,&
& '      +lambda_pw  : ',lambda_pw*1000._dp,ch10,&
& '      +lambda_paw : ',(lambda_v1-lambda_core-lambda_pw)*1000._dp,ch10
 call wrtout(std_out,msg,'COLL')
 write(msg,'(4a,es22.12,2a)') ch10,&
& ' Annihilation rate obtained from integration of e-p pairs momentum distribution:',ch10,&
& '   lambda=',lambda_v2*1000._dp,' ns-1',ch10
 call wrtout(ab_out,msg,'COLL')

!Deallocate remaining memory
 ABI_DEALLOCATE(my_gridtab)
 ABI_DEALLOCATE(rho_pw)
 ABI_DEALLOCATE(rho_moment_v1)
 ABI_DEALLOCATE(rho_moment_v2)
 ABI_DEALLOCATE(rho_moment_core)
 ABI_DEALLOCATE(rho_contrib)
 ABI_DEALLOCATE(rho_contrib_g)
 ABI_DEALLOCATE(rho_contrib_paw1)
 ABI_DEALLOCATE(rho_contrib_paw2)
 ABI_DEALLOCATE(rho_contrib_paw3)
 if (state_dependent) then
   call unset_mpi_enreg_fft(mpi_enreg_seq)
   call destroy_mpi_enreg(mpi_enreg_seq)
   ABI_DEALLOCATE(rhor_dop_el)
   if (dtset%usepaw==1) then
     call pawrhoij_free(pawrhoij_dop_el)
     ABI_DATATYPE_DEALLOCATE(pawrhoij_dop_el)
     if (mpi_enreg%my_natom<dtset%natom) then
       call pawrhoij_free(pawrhoij_all)
       call pawrhoij_free(pawrhoij_ep_all)
       ABI_DATATYPE_DEALLOCATE(pawrhoij_all)
       ABI_DATATYPE_DEALLOCATE(pawrhoij_ep_all)
     end if
   end if
 end if

 ABI_DEALLOCATE(gamma)

 if (dtset%usepaw==1.and.(.not.include_nhat_in_gamma)) then
   ABI_DEALLOCATE(rhor_)
   ABI_DEALLOCATE(rhor_ep_)
 end if

 if (dtset%usepaw==1) then
   ABI_DEALLOCATE(nphicor)
   ABI_DEALLOCATE(lmncmax)
   do itypat=1,dtset%ntypat
     if (allocated(phicor(itypat)%value)) then
       ABI_DEALLOCATE(phicor(itypat)%value)
     end if
     if (allocated(indlmncor(itypat)%value)) then
       ABI_DEALLOCATE(indlmncor(itypat)%value)
     end if
     if (allocated(radsumc(itypat)%value)) then
       ABI_DEALLOCATE(radsumc(itypat)%value)
     end if
     if (allocated(radsum1(itypat)%value)) then
       ABI_DEALLOCATE(radsum1(itypat)%value)
     end if
     if (allocated(radsum2(itypat)%value)) then
       ABI_DEALLOCATE(radsum2(itypat)%value)
     end if
     if (allocated(radsum3(itypat)%value)) then
       ABI_DEALLOCATE(radsum3(itypat)%value)
     end if
   end do
   do iatom=1,dtset%natom
     if (allocated(gammastate_c(iatom)%value)) then
       ABI_DEALLOCATE(gammastate_c(iatom)%value)
     end if
   end do
   ABI_DATATYPE_DEALLOCATE(phicor)
   ABI_DATATYPE_DEALLOCATE(indlmncor)
   ABI_DATATYPE_DEALLOCATE(radsumc)
   ABI_DATATYPE_DEALLOCATE(radsum1)
   ABI_DATATYPE_DEALLOCATE(radsum2)
   ABI_DATATYPE_DEALLOCATE(radsum3)
   ABI_DATATYPE_DEALLOCATE(gammastate_c)
 end if

 DBG_EXIT("COLL")

end subroutine posdoppler
!!***

!!****f* ABINIT/posratecore
!! NAME
!! posratecore
!!
!! FUNCTION
!! Calculate the annihilataion rate of a given core state
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | nspden=number of spin-density components
!!   | ntypat=number of atom types
!!   | paral_kgb=flag controlling (k,g,bands) parallelization
!!   | pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!   | usepaw=flag for PAW
!!  iatom= index of the current atom in posdoppler
!!  mesh_sizej= size of the radial mesh for the current atom in posdoppler
!!  mpi_enreg= information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  option= if 1, use gamma
!!          if 2, use IPM (gamma=1)
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(my_natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  rate= annihilation rate of a given core state needed for state dependent scheme for doppler broadening
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      posdoppler
!!
!! CHILDREN
!!      gammapositron,mkdenpos,nderiv_gen,pawdensities,pawxcsum,simp_gen
!!      xmpi_sum
!!
!! SOURCE

subroutine posratecore(dtset,electronpositron,iatom,my_natom,mesh_sizej,mpi_enreg,&
&                      option,pawang,pawrad,pawrhoij,pawrhoij_ep,&
&                      pawtab,rate,rhocorej)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iatom,my_natom,option,mesh_sizej
 real(dp),intent(out) :: rate
 type(dataset_type), intent(in) :: dtset
 type(electronpositron_type),pointer :: electronpositron
 type(MPI_type),intent(in) :: mpi_enreg
 type(pawang_type), intent(in) :: pawang
!arrays
 real(dp),intent(in) :: rhocorej(mesh_sizej)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*dtset%usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom*dtset%usepaw)
 type(pawrhoij_type),intent(in),target :: pawrhoij_ep(my_natom*dtset%usepaw)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer :: cplex,ierr,igamma,ii,ilm,ipt,ir
 integer :: itypat,iwarn,iwarnj,iwarnp,lm_size,lmn2_size,mesh_size
 integer :: ngr,ngrad,nspden_ep,opt_dens
 logical,parameter :: include_nhat_in_gamma=.false.
 real(dp),parameter :: delta=1.d-4
 real(dp) :: fact,fact2,intg
 real(dp) :: mpibuf,rdum,sqfpi
 character(len=500) :: msg
!arrays
 logical,allocatable :: lmselect(:),lmselect_ep(:),lmselect_dum(:)
 real(dp),parameter :: qphon(3)=(/zero,zero,zero/),lsign(2)=(/one,-one/)
 real(dp),allocatable :: d1gam(:,:),d2gam(:,:),ff(:),gam_(:,:,:),gamma(:,:),gammam(:,:),gg(:,:)
 real(dp),allocatable :: grhocore2(:),grhocor2_(:),grhoe2(:),grho2_(:)
 real(dp),allocatable :: nhat1(:,:,:),nhat1_ep(:,:,:)
 real(dp),allocatable :: rho_(:),rho_ep_(:),rho1(:,:,:),rho1_ep(:,:,:)
 real(dp),allocatable :: rhoarr1(:),rhoarr1_ep(:),rhoarr2(:)
 real(dp),allocatable :: rhocore(:),rhocor_(:)
 real(dp),allocatable :: rhosph(:),rhosph_ep(:),rhotot(:,:),rhotot_ep(:,:)
 real(dp),allocatable :: trho1(:,:,:),trho1_ep(:,:,:)
 real(dp),allocatable :: v1sum(:,:),v2sum(:,:,:)
 type(pawrhoij_type),pointer :: pawrhoij_ep_(:)

! *************************************************************************

 DBG_ENTER("COLL")

!Tests for developers
 if (.not.associated(electronpositron)) then
   msg='electronpositron variable must be associated!'
   MSG_BUG(msg)
 end if
!Constants
 fact=0.0
 cplex=1;nspden_ep=1
 ngrad=1;if (electronpositron%ixcpositron==3.or.electronpositron%ixcpositron==31) ngrad=2
 iwarn=0;iwarnj=0;iwarnp=1
 sqfpi=sqrt(four_pi)

!Compatibility tests
 if (electronpositron%particle==EP_NOTHING) then
   msg='Not valid for electronpositron%particle=NOTHING!'
   MSG_BUG(msg)
 end if

 if (dtset%usepaw==1) then
   if(dtset%pawxcdev==0.and.ngrad==2) then
     msg='GGA is not implemented for pawxcdev=0 (use dtset%pawxcdev/=0)!'
     MSG_BUG(msg)
   end if
 end if

!Select type(s) of enhancement factor
 if (electronpositron%ixcpositron==-1) igamma=0
 if (electronpositron%ixcpositron== 2) igamma=4
 if (electronpositron%ixcpositron==11.or.electronpositron%ixcpositron==31) igamma=3
 if (electronpositron%ixcpositron==1.or.electronpositron%ixcpositron==3) igamma=2
 if (option==2) igamma=0

 pawrhoij_ep_ => pawrhoij_ep

 rate=zero

 itypat=pawrhoij(iatom)%itypat
 lmn2_size=pawtab(itypat)%lmn2_size
 mesh_size=pawtab(itypat)%mesh_size
 lm_size=pawtab(itypat)%lcut_size**2
 cplex=1
 ngr=0;if (ngrad==2) ngr=mesh_size

!Allocations of "on-site" densities
 ABI_ALLOCATE(rho1 ,(cplex*mesh_size,lm_size,nspden_ep))
 ABI_ALLOCATE(trho1,(cplex*mesh_size,lm_size,nspden_ep))
 ABI_ALLOCATE(rho1_ep ,(cplex*mesh_size,lm_size,nspden_ep))
 ABI_ALLOCATE(trho1_ep,(cplex*mesh_size,lm_size,nspden_ep))
 ABI_ALLOCATE(lmselect,(lm_size))
 ABI_ALLOCATE(lmselect_ep,(lm_size))
 ABI_ALLOCATE(lmselect_dum,(lm_size))
 if (include_nhat_in_gamma) then
   ABI_ALLOCATE(nhat1,(cplex*mesh_size,lm_size,nspden_ep))
   ABI_ALLOCATE(nhat1_ep,(cplex*mesh_size,lm_size,nspden_ep))
 else
   ABI_ALLOCATE(nhat1,(0,0,0))
   ABI_ALLOCATE(nhat1_ep,(0,0,0))
 end if

!Compute "on-site" densities (n1, ntild1, nhat1) for electron and positron =====
 lmselect(:)=.true.
 opt_dens=1;if (include_nhat_in_gamma) opt_dens=0
 call pawdensities(rdum,cplex,iatom,lmselect,lmselect_dum,lm_size,nhat1,nspden_ep,1,&
& 0,opt_dens,-1,0,pawang,0,pawrad(itypat),pawrhoij(iatom),&
& pawtab(itypat),rho1,trho1)
 lmselect_ep(:)=.true.
 call pawdensities(rdum,cplex,iatom,lmselect_ep,lmselect_dum,lm_size,nhat1_ep,nspden_ep,1,&
& 0,opt_dens,-1,0,pawang,0,pawrad(itypat),pawrhoij_ep_(iatom),&
& pawtab(itypat),rho1_ep,trho1_ep)
!Compute contribution to annihilation rate

 ABI_ALLOCATE(rhocore,(mesh_size))

!First formalism: use densities on r,theta,phi
 if (dtset%pawxcdev==0) then

   ABI_ALLOCATE(gamma,(mesh_size,2))
   ABI_ALLOCATE(rhoarr1,(mesh_size))
   ABI_ALLOCATE(rhoarr1_ep,(mesh_size))

!  Loop on the angular part
   do ipt=1,pawang%angl_size
!    Build densities
     rhoarr1=zero;rhoarr1_ep=zero;rhocore=zero
     do ilm=1,lm_size
       if (lmselect(ilm)) rhoarr1(:)=rhoarr1(:)+rho1(:,ilm,1)*pawang%ylmr(ilm,ipt)
     end do
     do ilm=1,lm_size
       if (lmselect_ep(ilm)) rhoarr1_ep(:)=rhoarr1_ep(:)+rho1_ep(:,ilm,1)*pawang%ylmr(ilm,ipt)
     end do
     rhocore(:)=pawtab(itypat)%coredens(:)
!    Make the densities positive
     if (electronpositron%particle==EP_ELECTRON) then
       call mkdenpos(iwarnp,mesh_size,1,1,rhoarr1   ,dtset%xc_denpos)
       call mkdenpos(iwarn ,mesh_size,1,1,rhoarr1_ep,dtset%xc_denpos)
     else if (electronpositron%particle==EP_POSITRON) then
       call mkdenpos(iwarn ,mesh_size,1,1,rhoarr1   ,dtset%xc_denpos)
       call mkdenpos(iwarnp,mesh_size,1,1,rhoarr1_ep,dtset%xc_denpos)
     end if
!    Compute Gamma
     ABI_ALLOCATE(grhoe2,(ngr))
     ABI_ALLOCATE(grhocore2,(ngr))
     if (electronpositron%particle==EP_ELECTRON) then
       call gammapositron(gamma,grhocore2,grhoe2,igamma,ngr,mesh_size,&
&       rhocore,rhoarr1_ep,rhoarr1,1)
     else if (electronpositron%particle==EP_POSITRON) then
       call gammapositron(gamma,grhocore2,grhoe2,igamma,ngr,mesh_size,&
&       rhocore,rhoarr1,rhoarr1_ep,1)
     end if
     ABI_DEALLOCATE(grhoe2)
     ABI_DEALLOCATE(grhocore2)
!    Compute contribution to annihilation rates

     ABI_ALLOCATE(ff,(mesh_size))
     ff(1:mesh_size)=rhoarr1_ep(1:mesh_size)*rhocorej(1:mesh_size) &
&     *gamma(1:mesh_size,1)*pawrad(itypat)%rad(1:mesh_size)**2
     call simp_gen(intg,ff,pawrad(itypat))
     intg=intg*pawang%angwgth(ipt)*four_pi
     rate         =rate         +intg
     ABI_DEALLOCATE(ff)
   end do ! ipt
   ABI_DEALLOCATE(gamma)
   ABI_DEALLOCATE(rhoarr1)
   ABI_DEALLOCATE(rhoarr1_ep)

!Second formalism: use (l,m) moments for densities
 else if (dtset%pawxcdev/=0) then

!  Build densities
   ABI_ALLOCATE(gammam,(mesh_size,lm_size))
   ABI_ALLOCATE(rhotot,(mesh_size,lm_size))
   ABI_ALLOCATE(rhotot_ep,(mesh_size,lm_size))
   ABI_ALLOCATE(rhosph,(mesh_size))
   ABI_ALLOCATE(rhosph_ep,(mesh_size))

   rhotot   (:,:)=rho1   (:,:,1)
   rhotot_ep(:,:)=rho1_ep(:,:,1)
   rhocore(:)=pawtab(itypat)%coredens(:)
   rhosph   (:)=rhotot   (:,1)/sqfpi
   rhosph_ep(:)=rhotot_ep(:,1)/sqfpi
!  Make spherical densities positive
   if (electronpositron%particle==EP_ELECTRON) then
     call mkdenpos(iwarnp,mesh_size,1,1,rhosph   ,dtset%xc_denpos)
     call mkdenpos(iwarn ,mesh_size,1,1,rhosph_ep,dtset%xc_denpos)
   else if (electronpositron%particle==EP_POSITRON) then
     call mkdenpos(iwarn ,mesh_size,1,1,rhosph   ,dtset%xc_denpos)
     call mkdenpos(iwarnp,mesh_size,1,1,rhosph_ep,dtset%xc_denpos)
   end if

!  Need gradients of electronic densities for GGA
   ABI_ALLOCATE(grhoe2,(ngr))
   ABI_ALLOCATE(grhocore2,(ngr))
   if (ngr>0) then
     if (electronpositron%particle==EP_ELECTRON) then
       call nderiv_gen(grhoe2,rhosph_ep,pawrad(itypat))
     else if (electronpositron%particle==EP_POSITRON) then
       call nderiv_gen(grhoe2,rhosph,pawrad(itypat))
     end if
     grhoe2(:)=grhoe2(:)**2
     call nderiv_gen(grhocore2,rhocore,pawrad(itypat))
     grhocore2(:)=grhocore2(:)**2
   end if
!  Compute Gamma for (rho-,rho+),
!  (rho- +drho-,rho+), (rho- -drho-,rho+),
!  (rho-,rho+ +drho+), (rho-,rho+ -drho+),
!  (rho- +drho-,rho+ +drho+), (rho- -drho-,rho+ -drho+)
!  Do a seven steps loop
   ABI_ALLOCATE(gam_,(mesh_size,2,7))
   ABI_ALLOCATE(rho_,(mesh_size))
   ABI_ALLOCATE(rho_ep_,(mesh_size))
   ABI_ALLOCATE(rhocor_,(mesh_size))
   ABI_ALLOCATE(grho2_,(ngr))
   ABI_ALLOCATE(grhocor2_,(ngr))

   do ii=1,7
!    Apply delta to get perturbed densities
     rho_(:)=rhosph(:);rho_ep_(:)=rhosph_ep(:);rhocor_(:)=rhocore(:)
     if (ngr>0) grho2_(:)=grhoe2(:)
     if (ngr>0) grhocor2_(:)=grhocore2(:)
     if (ii==2.or.ii==4.or.ii==6) fact=(one+delta)
     if (ii==3.or.ii==5.or.ii==7) fact=(one-delta)
     fact2=fact**2
     if (ii==2.or.ii==3.or.ii==6.or.ii==7) then
       rho_(:)=fact*rho_(:)
       if (electronpositron%particle==EP_POSITRON) then
         if (ngr>0) grho2_(:)=fact2*grho2_(:)
         rhocor_(:)=fact*rhocor_(:)
         if (ngr>0) grhocor2_(:)=fact2*grhocor2_(:)
       end if
     end if

     if (ii==4.or.ii==5.or.ii==6.or.ii==7) then
       rho_ep_(:)=fact*rho_ep_(:)
       if (electronpositron%particle==EP_ELECTRON) then
         if (ngr>0) grho2_(:)=fact2*grho2_(:)
         rhocor_(:)=fact*rhocor_(:)
         if (ngr>0) grhocor2_(:)=fact2*grhocor2_(:)
       end if
     end if
!    Compute gamma for these perturbed densities
     if (electronpositron%particle==EP_ELECTRON) then
       call gammapositron(gam_(:,:,ii),grhocor2_,grho2_,igamma,ngr,mesh_size,rhocor_,rho_ep_,rho_,1)
     else if (electronpositron%particle==EP_POSITRON) then
       call gammapositron(gam_(:,:,ii),grhocor2_,grho2_,igamma,ngr,mesh_size,rhocor_,rho_,rho_ep_,1)
     end if

   end do ! end loop ii=1,7

   ABI_DEALLOCATE(rhocor_)
   ABI_DEALLOCATE(grho2_)
   ABI_DEALLOCATE(grhocor2_)
   ABI_DEALLOCATE(grhoe2)
   ABI_DEALLOCATE(grhocore2)
   rho_   (:)=rhosph   (:);if (electronpositron%particle==EP_POSITRON) rho_   (:)=rho_   (:)+rhocore(:)
   rho_ep_(:)=rhosph_ep(:);if (electronpositron%particle==EP_ELECTRON) rho_ep_(:)=rho_ep_(:)+rhocore(:)
!  Compute numerical first and second derivatives of Gamma
!  d1gam(1) = dgam/drho+ (particle=ELECTRON), dgam/drho- (particle=POSITRON)
!  d1gam(2) = dgam/drho- (particle=ELECTRON), dgam/drho+ (particle=POSITRON)
   ABI_ALLOCATE(d1gam,(mesh_size,2))
   d1gam(:,:)=zero
   do ir=1,mesh_size
     if (rho_     (ir)>tol14) d1gam(ir,1)=(gam_(ir,1,2)-gam_(ir,1,3))*half/(delta*rho_     (ir))
     if (rho_ep_  (ir)>tol14) d1gam(ir,2)=(gam_(ir,1,4)-gam_(ir,1,5))*half/(delta*rho_ep_  (ir))
   end do

!  d2gam(1) = d2gam/drho+_drho+ (particle=ELECTRON), dgam/drho-_drho- (particle=POSITRON)
!  d2gam(2) = d2gam/drho-_drho+ (particle=ELECTRON), dgam/drho+_drho- (particle=POSITRON)
!  d2gam(3) = d2gam/drho-_drho- (particle=ELECTRON), dgam/drho+_drho+ (particle=POSITRON)
   ABI_ALLOCATE(d2gam,(mesh_size,3))
   d2gam(:,:)=zero
   do ir=1,mesh_size
     if (rho_  (ir)>tol14) d2gam(ir,1)=(gam_(ir,1,2)+gam_(ir,1,3)-two*gam_(ir,1,1))/(delta*rho_  (ir))**2
     if (rho_ep_(ir)>tol14) then
       d2gam(ir,3)=(gam_(ir,1,4)+gam_(ir,1,5)-two*gam_(ir,1,1))/(delta*rho_ep_(ir))**2
       if (rho_(ir)>tol14) then
         d2gam(ir,2)=(gam_(ir,1,6)+gam_(ir,1,7)+two*gam_(ir,1,1) &
&         -gam_(ir,1,2)-gam_(ir,1,3)-gam_(ir,1,4)-gam_(ir,1,5)) &
&         *half/(delta*rho_(ir))/(delta*rho_ep_(ir))
       end if
     end if
   end do

   ABI_DEALLOCATE(rho_)
   ABI_DEALLOCATE(rho_ep_)
!  Compute useful sums of densities
   ABI_ALLOCATE(v1sum,(mesh_size,3))
   if ( dtset%pawxcdev>=2)  then
     ABI_ALLOCATE(v2sum,(mesh_size,lm_size,3))
   else
     ABI_ALLOCATE(v2sum,(0,0,0))
   end if
   rhotot(:,1)=sqfpi*rhosph(:);rhotot_ep(:,1)=sqfpi*rhosph_ep(:)
   call pawxcsum(1,1,1,lmselect,lmselect_ep,lm_size,mesh_size,3,dtset%pawxcdev,&
&   pawang,rhotot,rhotot_ep,v1sum,v2sum)
!  Compute final developpment of gamma moments
   gammam(:,:)=zero
   gammam(:,1)=gam_(:,1,1)*sqfpi
   gammam(:,1)=gammam(:,1)+(d2gam(:,2)*v1sum(:,2) &
&   +half*(d2gam(:,1)*v1sum(:,1)+d2gam(:,3)*v1sum(:,3)))/sqfpi
   do ilm=2,lm_size
     if (lmselect(ilm)) then
       gammam(:,ilm)=gammam(:,ilm)+d1gam(:,1)*rhotot(:,ilm)
     end if
     if (lmselect_ep(ilm)) then
       gammam(:,ilm)=gammam(:,ilm)+d1gam(:,2)*rhotot_ep(:,ilm)
     end if
   end do
   if (dtset%pawxcdev>1) then
     do ilm=2,lm_size
       gammam(:,ilm)=gammam(:,ilm)+d2gam(:,2)*v2sum(:,ilm,2) &
&       +half*(d2gam(:,1)*v2sum(:,ilm,1)+d2gam(:,3)*v2sum(:,ilm,3))
     end do
   end if

   ABI_DEALLOCATE(gam_)
   ABI_DEALLOCATE(d1gam)
   ABI_DEALLOCATE(d2gam)
   ABI_DEALLOCATE(v1sum)
   ABI_DEALLOCATE(v2sum)
!  Compute contribution to annihilation rate
   ABI_ALLOCATE(gg,(mesh_size,4))
   gg=zero
   ABI_ALLOCATE(rhoarr2,(mesh_size))
   do ilm=1,lm_size
     if (lmselect_ep(ilm)) gg(:,1)=gg(:,1)+rhotot_ep(:,ilm)*rhocorej(:)*gammam(:,ilm)
   end do
   ABI_DEALLOCATE(rhoarr2)
   gg(1:mesh_size,1)=gg(1:mesh_size,1)*pawrad(itypat)%rad(1:mesh_size)**2
   call simp_gen(intg,gg(:,1),pawrad(itypat))
   rate         =rate         +intg
   ABI_DEALLOCATE(gg)
   ABI_DEALLOCATE(gammam)
   ABI_DEALLOCATE(rhotot)
   ABI_DEALLOCATE(rhotot_ep)
   ABI_DEALLOCATE(rhosph)
   ABI_DEALLOCATE(rhosph_ep)

 end if ! dtset%pawxcdev
 ABI_DEALLOCATE(rhocore)

 ABI_DEALLOCATE(rho1)
 ABI_DEALLOCATE(trho1)
 ABI_DEALLOCATE(rho1_ep)
 ABI_DEALLOCATE(trho1_ep)
 ABI_DEALLOCATE(lmselect)
 ABI_DEALLOCATE(lmselect_ep)
 ABI_DEALLOCATE(lmselect_dum)
 ABI_DEALLOCATE(nhat1)
 ABI_DEALLOCATE(nhat1_ep)

!Reduction in case of distribution over atomic sites
 if (mpi_enreg%nproc_atom>1) then
   mpibuf=rate
   call xmpi_sum(mpibuf,mpi_enreg%comm_atom,ierr)
   rate=mpibuf
 end if

 DBG_EXIT("COLL")

end subroutine posratecore
!!***

end module m_positron
!!***
