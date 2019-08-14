!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_d2frnl
!! NAME
!!  m_d2frnl
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA, XG, GM, AR, MB, MT, AM)
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

module m_d2frnl

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_abicore
 use m_xmpi
 use m_errors
 use m_cgtools
 use m_nctk
 use m_hamiltonian
 use m_efmas_defs
 use m_wfk
 use m_dtset

 use m_time,     only : timab
 use m_geometry, only : metric, strconv
 use m_efmas,    only : check_degeneracies
 use m_io_tools, only : file_exists
 use m_hdr,      only : hdr_skip
 use m_pawang,   only : pawang_type
 use m_pawrad,   only : pawrad_type
 use m_pawtab,   only : pawtab_type,pawtab_get_lsize
 use m_pawfgrtab,only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_free
 use m_paw_ij,   only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify, paw_ij_reset_flags
 use m_pawrhoij, only : pawrhoij_type, pawrhoij_copy, pawrhoij_free, pawrhoij_gather, &
                        pawrhoij_nullify, pawrhoij_symrhoij
 use m_pawcprj,  only : pawcprj_type, pawcprj_alloc, pawcprj_get, pawcprj_copy, pawcprj_free
 use m_pawdij,   only : pawdijfr
 use m_paw_dfpt, only : pawgrnl
 use m_kg,       only : mkkin, mkkpg
 use m_mkffnl,   only : mkffnl
 use m_mpinfo,   only : proc_distrb_cycle
 use m_nonlop,   only : nonlop
 use m_paw_occupancies, only : pawaccrhoij

 implicit none

 private
!!***

 public :: d2frnl
!!***

contains
!!***

!!****f* ABINIT/d2frnl
!! NAME
!! d2frnl
!!
!! FUNCTION
!! Compute the frozen-wavefunction non-local contribution for response functions
!! (strain and/or phonon)
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=<G|Cnk>=Fourier coefficients of WF
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  dyfr_cplex=1 if dyfrnl is real, 2 if it is complex
!!  dyfr_nondiag=1 if dyfrnl is non diagonal with respect to atoms; 0 otherwise
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  gsqcut=Fourier cutoff on G^2 for "large sphere" of radius double that of the basis sphere
!!  has_allddk= True if all ddk file are present on disk
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  kg(3,mpw*mkmem)=work array for coordinates of G vectors in basis
!!   primitive translations
!!  mgfftf=maximum size of 1D FFTs for the fine FFT grid (PAW)
!!  mpi_enreg=information about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in unit cell
!!  nfftf= -PAW ONLY- number of FFT grid points for the fine grid
!!         (nfftf=nfft for norm-conserving potential runs)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/variables/vargs.htm#ngfft
!!  ngfftf(18)= -PAW ONLY- contain all needed information about 3D FFT for the fine grid
!!              (ngs_rbzfftf=ngfft for norm-conserving potential runs)
!!  npwarr(nkpt)=number of planewaves at each k point, and boundary
!!  ntypat=integer specification of atom type (1, 2, ...)
!!  occ(mband*nkpt*nsppol)=occupation numbers of bands (usually 2) at each k point
!!  rfphon=1   if non local contribution of dynamical matrix have to be computed
!!  rfstrs!=0  if non local contribution of elastic tensor have to be computed
!!  paw_ij(my_natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawbec= flag for the computation of Born Effective Charge within PAW ; set to 1 if yes
!!  pawpiezo= flag for the computation of piezoelectric tensor  within PAW ; set to 1 if yes
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=phase information related to structure factor
!!  ph1df(2,3*(2*mgfftf+1)*natom)=phase information related to structure factor on the fine FFT grid (PAW)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  symrec(3,3,nsym)=symmetries in reciprocal space (dimensionless)
!!  vtrial(nfftf,nspden)=total potential (Hartree+XC+loc)
!!  vxc(nfftf,nspden)=XC potential
!!  xred(3,natom)=reduced coordinates of atoms (dimensionless)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,9,mpsang*mpsang*useylm)= gradients of real spherical harmonics for each G and k point
!!
!! OUTPUT
!!  becfrnl(3,natom,3*pawbec)=NL frozen contribution to Born Effective Charges (PAW only)
!!                            (3,natom) = derivative wr to the displ. of one atom in one direction
!!                            (3)       = derivative wr to electric field in one direction
!!  piezofrnl(3,6*pawpiezo)=NL frozen contribution to piezoelectric tensor (PAW only)
!!  dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)=
!!         non-symmetrized non-local contribution to the dynamical matrix
!!         If NCPP, it depends on one atom
!!         If PAW,  it depends on two atoms
!!  eltfrnl(6+3*natom,6)=non-symmetrized non-local contribution to the
!!                    elastic tensor
!!
!! SIDE EFFECTS
!!  ===== if psps%usepaw==1
!!  pawfgrtab(my_natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!                          pawfgrtab(:)%gylmgr2 are deallocated here
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!    (gradients of rhoij for each atom with respect to atomic positions are computed here)
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      appdig,check_degeneracies,dotprod_g
!!      init_hamiltonian,metric,mkffnl
!!      mkkin,mkkpg,nonlop,paw_ij_free,paw_ij_init,paw_ij_nullify
!!      paw_ij_reset_flags,pawaccrhoij,pawcprj_alloc,pawcprj_free,pawdij2e1kb
!!      pawdijfr,pawfgrtab_free,pawfgrtab_init,pawgrnl,pawrhoij_free
!!      pawrhoij_gather,pawrhoij_nullify,pawtab_get_lsize,strconv,pawrhoij_symrhoij
!!      timab,wfk_close,wfk_open_read,wfk_read_bks,wrtout,xmpi_sum
!!
!! SOURCE

subroutine d2frnl(becfrnl,cg,dtfil,dtset,dyfrnl,dyfr_cplex,dyfr_nondiag,efmasdeg,efmasval,eigen,eltfrnl,&
&          gsqcut,has_allddk,indsym,kg,mgfftf,mpi_enreg,mpsang,my_natom,natom,nfftf,ngfft,ngfftf,npwarr,&
&          occ,paw_ij,pawang,pawbec,pawfgrtab,pawpiezo,pawrad,pawrhoij,pawtab,ph1d,ph1df,piezofrnl,psps,&
&          rprimd,rfphon,rfstrs,symrec,vtrial,vxc,xred,ylm,ylmgr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dyfr_cplex,dyfr_nondiag,mgfftf,mpsang,my_natom,natom
 integer,intent(in) :: nfftf,pawbec,pawpiezo,rfphon,rfstrs
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(in) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: indsym(4,dtset%nsym,natom),kg(3,dtset%mpw*dtset%mkmem)
 integer,intent(in) :: ngfft(18),ngfftf(18),npwarr(dtset%nkpt)
 integer,intent(in) :: symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: cg(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
 real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*natom)
 real(dp),intent(in) :: ph1df(2,3*(2*mgfftf+1)*natom),rprimd(3,3)
 real(dp),intent(in) :: vxc(nfftf,dtset%nspden),xred(3,natom)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,9,mpsang*mpsang*psps%useylm)
 real(dp),intent(in),target :: vtrial(nfftf,dtset%nspden)
 real(dp),intent(out) :: becfrnl(3,natom,3*pawbec),piezofrnl(6,3*pawpiezo)
 real(dp),intent(out) :: dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)
 real(dp),intent(out) :: eltfrnl(6+3*natom,6)
 logical,intent(inout):: has_allddk
 type(efmasdeg_type),allocatable,intent(out):: efmasdeg(:)
 type(efmasval_type),allocatable,intent(out):: efmasval(:,:)
 type(paw_ij_type),intent(in) :: paw_ij(my_natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat)
 type(pawrhoij_type),intent(inout),target :: pawrhoij(my_natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat)

!Local variables-------------------------------
!scalars
 integer,parameter :: formeig1=1,usecprj=0
 integer :: bandmin,bandmax,bdtot_index,bufdim
 integer :: choice_bec2,choice_bec54,choice_efmas,choice_phon,choice_strs,choice_piez3,choice_piez55
 integer :: cplex,cplx,cpopt,cpopt_bec,ddkcase,deg_dim
 integer :: dimffnl,dimffnl_str,dimnhat,ia,iatom,iashift,iband,jband,ibg,icg,icplx,ideg,ider,idir
 integer :: ider_str,idir_ffnl,idir_str,ielt,ieltx,ierr,ii,ikg,ikpt,ilm,ipw,iq,iq0
 integer :: ispinor,isppol,istwf_k,isub,itypat,jj,jsub,klmn,master,me,mu
 integer :: my_comm_atom,n1,n2,n3,nband_k,ncpgr,nfftot,ngrhoij,nkpg,nnlout_bec1,nnlout_bec2,nnlout_efmas
 integer :: nnlout_piez1,nnlout_piez2,nnlout_phon,nnlout_strs,npw_,npw_k,nsp,nsploop,nu
 integer :: optgr,optgr2,option,option_rhoij,optstr,optstr2,paw_opt,paw_opt_1,paw_opt_3,paw_opt_efmas
 integer :: shift_rhoij,signs,signs_field,spaceworld,sz2,sz3,tim_nonlop
 real(dp) :: arg,eig_k,enl,enlk,occ_k,ucvol,wtk_k
 logical :: has_ddk_file,need_becfr,need_efmas,need_piezofr,paral_atom,t_test,usetimerev
 character(len=500) :: msg
 type(gs_hamiltonian_type) :: gs_ham
!arrays
 integer :: ik_ddk(3),ddkfil(3)
 integer,allocatable :: dimlmn(:),kg_k(:,:),l_size_atm(:)
 integer,pointer :: my_atmtab(:)
 real(dp) :: dotprod(2),dummy(0),gmet(3,3),gprimd(3,3),grhoij(3),kpoint(3),nonlop_dum(1,1)
 real(dp) :: rmet(3,3),tsec(2)
 real(dp),allocatable :: becfrnl_tmp(:,:,:),becfrnlk(:,:,:),becij(:,:,:,:,:),cg_left(:,:)
 real(dp),allocatable :: cwavef(:,:),ddk(:,:),ddkinpw(:,:,:),dyfrnlk(:,:)
 real(dp),allocatable :: elt_work(:,:),eltfrnlk(:,:),enlout_bec1(:),enlout_bec2(:),enlout_efmas(:)
 real(dp),allocatable :: enlout_piez1(:),enlout_piez2(:),enlout_phon(:),enlout_strs(:)
 real(dp),allocatable :: gh2c(:,:),gs2c(:,:)
 real(dp),allocatable :: kpg_k(:,:),mpibuf(:),nhat_dum(:,:),piezofrnlk(:,:),ph3d(:,:,:)
 real(dp),allocatable :: svectout(:,:),ylm_k(:,:),ylmgr_k(:,:,:)
 real(dp),allocatable,target :: ffnl(:,:,:,:),ffnl_str(:,:,:,:,:)
 character(len=fnlen) :: fiwfddk(3)
 type(paw_ij_type),allocatable :: paw_ij_tmp(:)
 type(pawcprj_type),allocatable,target :: cwaveprj(:,:)
 type(pawfgrtab_type),allocatable :: pawfgrtab_tmp(:)
 type(pawrhoij_type),pointer :: pawrhoij_tot(:)
 type(wfk_t) :: ddkfiles(3)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(159,1,tsec)

 write(msg,'(3a)')ch10,' ==> Calculation of the frozen part of the second order derivatives, this can take some time...',ch10
 call wrtout(std_out,msg,'COLL')

!Set up parallelism
 spaceworld=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt
 master=0
 paral_atom=(my_natom/=natom)
 my_comm_atom=mpi_enreg%comm_atom
 my_atmtab=>mpi_enreg%my_atmtab

!Compute gmet, gprimd and ucvol from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!If needed, check for ddk files (used for effective charges)
 if (pawbec==1.or.pawpiezo==1) then
   ddkfil(:)=0
   do ii=1,3
     ddkcase=ii+natom*3
     call appdig(ddkcase,dtfil%fnamewffddk,fiwfddk(ii))
     t_test = file_exists(fiwfddk(ii))
     ! Trick needed to run Abinit test suite in netcdf mode.
     if (.not. t_test .and. file_exists(nctk_ncify(fiwfddk(ii)))) then
       t_test = .True.; fiwfddk(ii) = nctk_ncify(fiwfddk(ii))
       write(msg,"(3a)")"- File: ",trim(fiwfddk(ii))," does not exist but found netcdf file with similar name."
       call wrtout(std_out,msg,'COLL')
     end if
     if (t_test) ddkfil(ii)=20+ii ! Note the use of unit numbers 21, 22 and 23
   end do
   has_ddk_file=(any(ddkfil(:)>0))
   has_allddk  =(all(ddkfil(:)>0))
 else
   has_ddk_file=.FALSE.
   has_allddk  =.FALSE.
 end if

 if(pawbec==1.or.pawpiezo==1.and.has_ddk_file) then
   if(.not.has_allddk) then
     write(msg,'(5a)')ch10,&
&     ' WARNING: All ddk perturbations are needed to compute',ch10,&
&     ' the frozen part of Born effective charges and/or piezoelectric tensor.',ch10
     call wrtout(std_out,msg,'COLL')
   else
     write(msg,'(5a)')ch10,&
&     ' All ddk perturbations are available.',ch10,&
&     ' The frozen part of Born effective charges and/or piezoelectric tensor will be computed',ch10
     call wrtout(std_out,msg,'COLL')
   end if
 end if

 need_becfr=(pawbec==1.and.has_ddk_file)
 need_piezofr=(pawpiezo==1.and.has_ddk_file)

!Initialization of frozen non local array
 if(rfphon==1) then
   dyfrnl(:,:,:,:,:)=zero
   ABI_ALLOCATE(dyfrnlk,(6,natom))
 end if
 if(rfstrs/=0)then
   eltfrnl(:,:)=zero;enl=zero
   ABI_ALLOCATE(eltfrnlk,(6+3*natom,6))
 end if
 if (need_becfr) then
   becfrnl(:,:,:)=zero
   ABI_ALLOCATE(becfrnlk,(3,natom,3))
 end if
 if (need_piezofr) then
   piezofrnl(:,:)=zero
   ABI_ALLOCATE(piezofrnlk,(6,3))
 end if
 need_efmas=dtset%efmas>0
 if(need_efmas.and.(rfphon==1.or.rfstrs/=0.or.need_becfr.or.need_piezofr)) then
   write(msg,'(5a)')ch10,&
&   ' ERROR: Efmas calculation is incompatible with phonons, elastic tensor, Born effective charges,',ch10,&
&   ' and piezoelectric tensor calculations. Please revise your input.',ch10
   MSG_ERROR(msg)
 end if

!Common initialization
 bdtot_index=0;ibg=0;icg=0
 nsploop=dtset%nsppol;if (dtset%nspden==4) nsploop=4
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 nfftot=ngfftf(1)*ngfftf(2)*ngfftf(3)

!Common data for "nonlop" routine
 tim_nonlop=6
 signs=1 ; signs_field = 2 ; eig_k=zero ; idir=0
! ffnl are in cartesian coordinates for EFMAS (idir==4),
! in contrast to reduced coordinates for the other responses (idir==0).
 idir_ffnl=0 ; if(need_efmas) idir_ffnl=4
 choice_phon=0;choice_strs=0
 if(rfphon==1)then
   shift_rhoij=0
   choice_phon=4
   nnlout_phon=max(1,6*natom)
   ABI_ALLOCATE(enlout_phon,(nnlout_phon))
 end if
 if(rfstrs/=0)then
   shift_rhoij=6
   choice_strs=6
   nnlout_strs=6*(3*natom+6)
   ABI_ALLOCATE(enlout_strs,(nnlout_strs))
 end if
 if (psps%usepaw==0) then
   paw_opt=0 ; cpopt=-1
 else
   paw_opt=2 ; cpopt=1+2*usecprj
 end if
 if(need_piezofr)then
   choice_piez3  =  3
   choice_piez55 = 55
   nnlout_piez1  =  6
   nnlout_piez2  = 36
   paw_opt_1     = 1
   paw_opt_3     = 3
   ABI_ALLOCATE(enlout_piez1,(nnlout_piez1))
   ABI_ALLOCATE(enlout_piez2,(nnlout_piez2))
 end if
 if (need_becfr) then
   choice_bec2=2 ; choice_bec54=54
   nnlout_bec1=max(1,3*natom) ; nnlout_bec2=max(1,18*natom);
   paw_opt_1=1 ; paw_opt_3=3 ; cpopt_bec=-1
   ABI_ALLOCATE(enlout_bec1,(nnlout_bec1))
   ABI_ALLOCATE(enlout_bec2,(nnlout_bec2))
 else
   choice_bec2=0 ; choice_bec54=0
   nnlout_bec1=0
 end if
 if(need_efmas) then
   ABI_MALLOC(enlout_efmas,(0))
   ABI_DATATYPE_ALLOCATE(efmasdeg,(dtset%nkpt))
   ABI_DATATYPE_ALLOCATE(efmasval,(dtset%mband,dtset%nkpt))
 end if

!Initialize Hamiltonian (k-independent terms)
 call init_hamiltonian(gs_ham,psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,natom,&
& dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,&
& paw_ij=paw_ij,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab,&
& usecprj=usecprj,ph1d=ph1d,nucdipmom=dtset%nucdipmom)

!===== PAW specific section
 if (psps%usepaw==1) then

!  Define several sizes & flags
   ncpgr=0;ngrhoij=0
   if(rfphon==1)then
     ncpgr=3;ngrhoij=3
   end if
   if(rfphon==1.or.need_becfr)then
     ncpgr=6;ngrhoij=6
   end if
   if(rfstrs/=0.and.rfphon==1)then
     ncpgr=9;ngrhoij=9
   end if
   if(rfstrs/=0.or.need_piezofr)then
     ncpgr=9;ngrhoij=9
   end if

!  If PAW and Born Eff. Charges, one has to compute some additional data:
!  For each atom and for electric field direction k:
!  becij(k)=<Phi_i|r_k-R_k|Phi_j>-<tPhi_i|r_k-R_k|tPhi_j> + sij.R_k
   if (need_becfr.or.need_piezofr) then
     ABI_ALLOCATE(becij,(gs_ham%dimekb1,gs_ham%dimekb2,dtset%nspinor**2,1,3))
     becij=zero
     ABI_DATATYPE_ALLOCATE(paw_ij_tmp,(my_natom))
     ABI_DATATYPE_ALLOCATE(pawfgrtab_tmp,(my_natom))
     call paw_ij_nullify(paw_ij_tmp)
     cplex=1;nsp=1 ! Force nsppol/nspden to 1 because Dij^(1) due to electric field is spin-independent
     call paw_ij_init(paw_ij_tmp,cplex,dtset%nspinor,nsp,nsp,dtset%pawspnorb,natom,psps%ntypat,&
&     dtset%typat,pawtab,has_dijfr=1,comm_atom=my_comm_atom,mpi_atmtab=my_atmtab )
     call pawtab_get_lsize(pawtab,l_size_atm,my_natom,dtset%typat,mpi_atmtab=my_atmtab)
     call pawfgrtab_init(pawfgrtab_tmp,1,l_size_atm,dtset%nspden,dtset%typat,&
&     mpi_atmtab=my_atmtab,comm_atom=my_comm_atom)
     ABI_DEALLOCATE(l_size_atm)
     do ii=1,3 ! Loop over direction of electric field
       call paw_ij_reset_flags(paw_ij_tmp,all=.True.)
       call pawdijfr(gprimd,ii,natom+2,my_natom,natom,nfftf,ngfftf,nsp,nsp,psps%ntypat,&
&       0,paw_ij_tmp,pawang,pawfgrtab_tmp,pawrad,pawtab,cplex,&
&       (/zero,zero,zero/),rprimd,ucvol,vtrial,vtrial,vxc,xred,&
&       comm_atom=my_comm_atom, mpi_atmtab=my_atmtab ) ! vtrial not used here
       do isppol=1,dtset%nspinor**2
         call pawdij2e1kb(paw_ij_tmp(:),nsp,my_comm_atom,e1kbfr=becij(:,:,:,:,ii),mpi_atmtab=my_atmtab)
       end do
     end do
     call paw_ij_free(paw_ij_tmp)
     call pawfgrtab_free(pawfgrtab_tmp)
     ABI_DATATYPE_DEALLOCATE(paw_ij_tmp)
     ABI_DATATYPE_DEALLOCATE(pawfgrtab_tmp)
   end if

!  PAW occupancies: need to communicate when paral atom is activated
   if (paral_atom) then
     ABI_DATATYPE_ALLOCATE(pawrhoij_tot,(natom))
     call pawrhoij_nullify(pawrhoij_tot)
     call pawrhoij_gather(pawrhoij,pawrhoij_tot,-1,my_comm_atom)
   else
     pawrhoij_tot => pawrhoij
   end if

!  Projected WF (cprj) and PAW occupancies (& gradients)
   ABI_DATATYPE_ALLOCATE(cwaveprj,(natom,dtset%nspinor))
   call pawcprj_alloc(cwaveprj,ncpgr,gs_ham%dimcprj)
   do iatom=1,natom
     sz2=pawrhoij_tot(iatom)%cplex_rhoij*pawrhoij_tot(iatom)%qphase*pawrhoij_tot(iatom)%lmn2_size
     sz3=pawrhoij_tot(iatom)%nspden
     ABI_ALLOCATE(pawrhoij_tot(iatom)%grhoij,(ngrhoij,sz2,sz3))
     pawrhoij_tot(iatom)%ngrhoij=ngrhoij
     pawrhoij_tot(iatom)%grhoij=zero
   end do
   usetimerev=(dtset%kptopt>0.and.dtset%kptopt<3)

 else
   ABI_DATATYPE_ALLOCATE(cwaveprj,(0,0))
 end if !PAW

!If needed, manage ddk files
!Open ddk WF file(s)
 if (need_becfr.or.need_piezofr) then
   do ii=1,3 ! Loop over elect. field directions
     if (ddkfil(ii)/=0) then
       write(msg, '(a,a)') '-open ddk wf file :',trim(fiwfddk(ii))
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out,msg,'COLL')
       call wfk_open_read(ddkfiles(ii),fiwfddk(ii),formeig1,dtset%iomode,ddkfil(ii),spaceworld)
     end if
   end do
 end if

!LOOP OVER SPINS
 do isppol=1,dtset%nsppol

!  Continue to initialize the Hamiltonian (PAW DIJ coefficients)
   call gs_ham%load_spin(isppol,with_nonlocal=.true.)

!  Rewind (k+G) data if needed
   ikg=0

!  Loop over k points
   do ikpt=1,dtset%nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
     istwf_k=dtset%istwfk(ikpt)
     npw_k=npwarr(ikpt)
     wtk_k=dtset%wtk(ikpt)
     kpoint(:)=dtset%kptns(:,ikpt)

!    Skip this k-point if not the proper processor
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
       bdtot_index=bdtot_index+nband_k
       cycle
     end if

!    If needed, manage ddk files
     if (need_becfr.or.need_piezofr) then
       do ii=1,3 ! Loop over elect. field directions
         if (ddkfil(ii)/=0)then
!        Number of k points to skip in the full set of k pointsp
           ik_ddk(ii) = ddkfiles(ii)%findk(kpoint)
           ABI_CHECK(ik_ddk(ii) /= -1, "Cannot find k-point in DDK")
           npw_ = ddkfiles(ii)%hdr%npwarr(ik_ddk(ii))
           if (npw_/=npw_k) then
             write(unit=msg,fmt='(a,i3,a,i5,a,i3,a,a,i5,a,a,i5)')&
&             'For isppol = ',isppol,', ikpt = ',ikpt,' and idir = ',ii,ch10,&
&             'the number of plane waves in the ddk file is equal to', npw_,ch10,&
&             'while it should be ',npw_k
             MSG_ERROR(msg)
           end if

         end if
       end do
     end if

     ABI_ALLOCATE(cwavef,(2,npw_k*dtset%nspinor))
     if (need_becfr.or.need_piezofr) then
       ABI_ALLOCATE(svectout,(2,npw_k*dtset%nspinor))
     end if
     if (need_efmas) then
       ABI_MALLOC(cg_left,(2,npw_k*dtset%nspinor))
       ABI_MALLOC(gh2c,(2,npw_k*dtset%nspinor))
       ABI_MALLOC(gs2c,(2,npw_k*dtset%nspinor))
     end if

     ABI_ALLOCATE(ylm_k,(npw_k,mpsang*mpsang*psps%useylm))
     if(rfstrs/=0.or.need_becfr.or.need_piezofr.or.need_efmas)then
       ABI_ALLOCATE(ylmgr_k,(npw_k,9,mpsang*mpsang*psps%useylm))
     else
       ABI_ALLOCATE(ylmgr_k,(0,0,0))
     end if

     ABI_ALLOCATE(kg_k,(3,npw_k))
     kg_k(:,:) = 0
!$OMP PARALLEL DO
     do ipw=1,npw_k
       kg_k(1,ipw)=kg(1,ipw+ikg)
       kg_k(2,ipw)=kg(2,ipw+ikg)
       kg_k(3,ipw)=kg(3,ipw+ikg)
     end do
     if (psps%useylm==1) then
!SOMP PARALLEL DO COLLAPSE(2)
       do ilm=1,mpsang*mpsang
         do ipw=1,npw_k
           ylm_k(ipw,ilm)=ylm(ipw+ikg,ilm)
         end do
       end do
       if(rfstrs/=0.or.need_becfr.or.need_piezofr.or.need_efmas)then
!SOMP PARALLEL DO COLLAPSE(3)
         do ilm=1,mpsang*mpsang
           do ii=1,9
             do ipw=1,npw_k
               ylmgr_k(ipw,ii,ilm)=ylmgr(ipw+ikg,ii,ilm)
             end do
           end do
         end do
       end if
     end if

     cplex=2;if (istwf_k>1) cplex=1

!    Compute (k+G) vectors (only if useylm=1)
     nkpg=0
     if (rfstrs/=0.or.need_efmas.or.pawpiezo==1) nkpg=3*dtset%nloalg(3)
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     if (nkpg>0) then
       call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
     end if

     !EFMAS: Compute second order derivatives w/r to k for all direction for this k-point.
     if (need_efmas) then
       ABI_ALLOCATE(ddkinpw,(npw_k,3,3))
       do mu=1,3
         do nu=1,3
!           call d2kpg(ddkinpw(:,mu,nu),dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,mu,nu,kg_k,kpoint,npw_k)
           call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,ddkinpw(:,mu,nu),kpoint,npw_k,mu,nu)
         end do
       end do
     end if

!    Compute nonlocal form factors ffnl at all (k+G):
     ider=0;dimffnl=1;
     if(need_becfr) then
       ider=1;dimffnl=4
     end if
     if(rfstrs/=0.or.need_piezofr.or.need_efmas)then
       ider=2;dimffnl=3+7*psps%useylm
     end if
     ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,psps%ntypat))
     call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&     gmet,gprimd,ider,idir_ffnl,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&     psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,&
&     psps%ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_k)

!    For piezoelectric tensor need additional ffnl derivatives
     if(need_piezofr)then
       ider_str=1 ; dimffnl_str=2
       ABI_ALLOCATE(ffnl_str,(npw_k,dimffnl_str,psps%lmnmax,psps%ntypat,6))
       do mu=1,6 !loop over strain
         idir_str=-mu
         call mkffnl(psps%dimekb,dimffnl_str,psps%ekb,ffnl_str(:,:,:,:,mu),&
&         psps%ffspl,gmet,gprimd,ider_str,idir_str,psps%indlmn,kg_k,kpg_k,&
&         kpoint,psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,&
&         psps%ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_k)
       end do
     end if

!    Load k-dependent part in the Hamiltonian datastructure
     ABI_ALLOCATE(ph3d,(2,npw_k,gs_ham%matblk))
     call gs_ham%load_k(kpt_k=kpoint,npw_k=npw_k,istwf_k=istwf_k,&
&     kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnl,ph3d_k=ph3d,compute_ph3d=.true.)

!    Initialize contributions from current k point
     if(rfphon==1) dyfrnlk(:,:)=zero
     if(rfstrs/=0)then
       enlk=zero;eltfrnlk(:,:)=zero
     end if
     if (need_becfr) becfrnlk(:,:,:)=zero
     if (need_piezofr) piezofrnlk(:,:)=zero
     if(need_efmas) then
       call check_degeneracies(efmasdeg(ikpt),dtset%efmas_bands(:,ikpt),nband_k,eigen(bdtot_index+1:bdtot_index+nband_k), &
&       dtset%efmas_deg_tol)
       do ideg=1,efmasdeg(ikpt)%ndegs
         if( efmasdeg(ikpt)%deg_range(1) <= ideg .and. ideg <= efmasdeg(ikpt)%deg_range(2) ) then
           deg_dim=efmasdeg(ikpt)%degs_bounds(2,ideg) - efmasdeg(ikpt)%degs_bounds(1,ideg) + 1
           ABI_MALLOC(efmasval(ideg,ikpt)%ch2c,(3,3,deg_dim,deg_dim))
           ABI_MALLOC(efmasval(ideg,ikpt)%eig2_diag,(3,3,deg_dim,deg_dim))
           efmasval(ideg,ikpt)%ch2c=zero
           efmasval(ideg,ikpt)%eig2_diag=zero
         else
           ABI_MALLOC(efmasval(ideg,ikpt)%ch2c,(0,0,0,0))
           ABI_MALLOC(efmasval(ideg,ikpt)%eig2_diag,(0,0,0,0))
         end if
       end do
     end if

!    Loop over bands
     do iband=1,nband_k

       if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband,iband,isppol,me)) then
         cycle
       end if

       occ_k=occ(iband+bdtot_index)
       cwavef(:,1:npw_k*dtset%nspinor) = cg(:,1+(iband-1)*npw_k*dtset%nspinor+icg:iband*npw_k*dtset%nspinor+icg)

!      Compute non-local contributions from n,k
       if (psps%usepaw==1) eig_k=eigen(iband+bdtot_index)

!      === Dynamical matrix
       if(rfphon==1) then
         call nonlop(choice_phon,cpopt,cwaveprj,enlout_phon,gs_ham,idir,(/eig_k/),mpi_enreg,1,&
&         nnlout_phon,paw_opt,signs,nonlop_dum,tim_nonlop,cwavef,cwavef)
!        Accumulate non-local contributions from n,k
         dyfrnlk(:,:)=dyfrnlk(:,:)+occ_k*reshape(enlout_phon(:),(/6,natom/))
       end if

!      === Elastic tensor
       if(rfstrs/=0) then
         call nonlop(choice_strs,cpopt,cwaveprj,enlout_strs,gs_ham,idir,(/eig_k/),mpi_enreg,1,&
&         nnlout_strs,paw_opt,signs,nonlop_dum,tim_nonlop,cwavef,cwavef)
!        Accumulate non-local contribut ions from n,k
         eltfrnlk(:,:)=eltfrnlk(:,:)+occ_k*reshape(enlout_strs(:),(/3*natom+6,6/))
       end if !endo if strs

!      PAW: accumulate gradients of rhoij
       !EFMAS: Bug with efmas currently; to be looked into...
       if (psps%usepaw==1.and.(.not.need_efmas)) then
         call pawaccrhoij(gs_ham%atindx,cplex,cwaveprj,cwaveprj,0,isppol,natom,&
&         natom,dtset%nspinor,occ_k,3,pawrhoij_tot,usetimerev,wtk_k)
       end if

!      PAW: Compute frozen contribution to piezo electric tensor
       if (need_piezofr) then
         do ii=1,3 ! Loop over elect. field directions
           call nonlop(choice_piez3,cpopt,cwaveprj,enlout_piez1,gs_ham,0,(/zero/),mpi_enreg,1,&
&           nnlout_piez1,paw_opt_1,signs,nonlop_dum,tim_nonlop,cwavef,cwavef,enl=becij(:,:,:,:,ii))
           piezofrnlk(:,ii)=piezofrnlk(:,ii)+occ_k*enlout_piez1(:)
         end do !end do ii
       end if

!      PAW: Compute frozen contribution to Born Effective Charges
       if (need_becfr) then
         do ii=1,3 ! Loop over elect. field directions
           call nonlop(choice_bec2,cpopt,cwaveprj,enlout_bec1,gs_ham,0,(/zero/),mpi_enreg,1,&
&           nnlout_bec1,paw_opt_1,signs,nonlop_dum,tim_nonlop,cwavef,cwavef,enl=becij(:,:,:,:,ii))
           becfrnlk(:,:,ii)=becfrnlk(:,:,ii)+occ_k*reshape(enlout_bec1(:),(/3,natom/))
         end do !end do ii
       end if

       if (need_becfr.or.need_piezofr) then
         do ii=1,3 ! Loop over elect. field directions
!          Not able to compute if ipert=(Elect. field) and no ddk WF file
           if (ddkfil(ii)==0) cycle
!            Read ddk wave function
           ABI_ALLOCATE(ddk,(2,npw_k*dtset%nspinor))
           if (ddkfil(ii)/=0) then
             call ddkfiles(ii)%read_bks(iband, ik_ddk(ii), isppol, xmpio_single, cg_bks=ddk)
!            Multiply ddk by +i
             do jj=1,npw_k*dtset%nspinor
               arg=ddk(1,jj)
               ddk(1,jj)=-ddk(2,jj);ddk(2,jj)=arg
             end do
           else
             ddk=zero
           end if

           if(need_becfr)then
             do iatom=1,natom !Loop over atom
               ia=gs_ham%atindx(iatom)
               do mu=1,3 !loop over atom direction
                 call nonlop(choice_bec2,cpopt_bec,cwaveprj,enlout_bec1,gs_ham,mu,(/zero/),&
&                 mpi_enreg,1,nnlout_bec1,paw_opt_3,signs_field,svectout,tim_nonlop,&
&                 cwavef,cwavef,iatom_only=iatom)
                 call dotprod_g(dotprod(1),dotprod(2),istwf_k,npw_k*dtset%nspinor,2,svectout,ddk,&
&                 mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
                 becfrnlk(mu,ia,ii)=becfrnlk(mu,ia,ii)+occ_k*dotprod(1)
               end do
             end do
           end if

           if(need_piezofr)then
             do mu=1,6 !loop over strain
               call gs_ham%load_k(ffnl_k=ffnl_str(:,:,:,:,mu))
               call nonlop(choice_piez3,cpopt,cwaveprj,enlout_piez1,gs_ham,mu,(/zero/),mpi_enreg,1,&
&               nnlout_piez1,paw_opt_3,signs_field,svectout,tim_nonlop,cwavef,svectout)
               call dotprod_g(dotprod(1),dotprod(2),istwf_k,npw_k*dtset%nspinor,2,svectout,ddk,&
&               mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
               piezofrnlk(mu,ii)=piezofrnlk(mu,ii)+occ_k*dotprod(1)
             end do
             call gs_ham%load_k(ffnl_k=ffnl)
           end if

           ABI_DEALLOCATE(ddk)
         end do ! End loop ddk file
       end if

       if(need_piezofr)then
         enlout_piez2 = zero
         call nonlop(choice_piez55,cpopt,cwaveprj,enlout_piez2,gs_ham,0,(/zero/),mpi_enreg,1,&
&         nnlout_piez2,paw_opt_3,signs,nonlop_dum,tim_nonlop,cwavef,cwavef)
!         Multiply enlout by +i
         iashift = 1
         do mu=1,6     ! strain
           do nu=1,3   ! k
             piezofrnlk(mu,nu)=piezofrnlk(mu,nu)-occ_k*(enlout_piez2(iashift+1)) ! Real part
!            piezofrnlk(mu,nu)=piezofrnlk(mu,nu)+occ_k*(enlout_piez2(iashift  ))! Imaginary part
             iashift = iashift + 2
           end do
         end do
       end if

       if(need_becfr)then
         call nonlop(choice_bec54,cpopt,cwaveprj,enlout_bec2,gs_ham,0,(/zero/),mpi_enreg,1,&
&         nnlout_bec2,paw_opt_3,signs,nonlop_dum,tim_nonlop,cwavef,cwavef)
!        Multiply enlout by +i
         iashift = 1
         do iatom=1,natom ! atm
           do mu=1,3     ! atm pos.
             do nu=1,3   ! k
               becfrnlk(mu,iatom,nu)=becfrnlk(mu,iatom,nu)-occ_k*(enlout_bec2(iashift+1)) ! Real part
!               becfrnlk(mu,iatom,nu)=becfrnlk(mu,iatom,nu)+occ_k*(enlout_bec2(iashift  ))! Imaginary part
               iashift = iashift + 2
             end do
           end do
         end do
       end if

       if(need_efmas) then
         bandmin=efmasdeg(ikpt)%degs_bounds(1, efmasdeg(ikpt)%deg_range(1) )
         bandmax=efmasdeg(ikpt)%degs_bounds(2, efmasdeg(ikpt)%deg_range(2) )
         if ( iband>=bandmin .and. iband<=bandmax ) then
           choice_efmas=8; signs=2
           cpopt=-1  !To prevent re-use of stored dgxdt, which are not for all direction required for EFMAS.
           paw_opt_efmas=0; if(psps%usepaw/=0) paw_opt_efmas=4 !To get both gh2c and gs2c
           nnlout_efmas=0; tim_nonlop=0 ! No tim_nonlop for efmas, currently.
           do mu=1,3
             do nu=1,3
               idir=3*(mu-1)+nu !xx=1, xy=2, xz=3, yx=4, yy=5, yz=6, zx=7, zy=8, zz=9, (xyz,xyz)=(mu,nu)
               gh2c=zero; gs2c=zero
               call nonlop(choice_efmas,cpopt,cwaveprj,enlout_efmas,gs_ham,idir,(/eig_k/),mpi_enreg,&
               1,nnlout_efmas,paw_opt_efmas,signs,gs2c,tim_nonlop,cwavef,gh2c)
               do ispinor=1,dtset%nspinor
                 ii = 1+(ispinor-1)*npw_k
                 do icplx=1,2
                   gh2c(icplx,ii:ispinor*npw_k) = gh2c(icplx,ii:ispinor*npw_k) +  &
&                   ddkinpw(1:npw_k,mu,nu)*cwavef(icplx,ii:ispinor*npw_k)
                 end do
               end do
               gh2c = gh2c - eig_k*gs2c
               ideg = efmasdeg(ikpt)%ideg(iband)
               do jband=efmasdeg(ikpt)%degs_bounds(1,ideg),efmasdeg(ikpt)%degs_bounds(2,ideg)
                 cg_left(:,1:npw_k*dtset%nspinor) = cg(:,1+(jband-1)*npw_k*dtset%nspinor+icg:jband*npw_k*dtset%nspinor+icg)
                 dotprod=0
                 call dotprod_g(dotprod(1),dotprod(2),istwf_k,npw_k*dtset%nspinor,2,cg_left,gh2c,mpi_enreg%me_g0,&
&                 mpi_enreg%comm_spinorfft)
                 isub = iband-efmasdeg(ikpt)%degs_bounds(1,ideg)+1
                 jsub = jband-efmasdeg(ikpt)%degs_bounds(1,ideg)+1
                 efmasval(ideg,ikpt)%ch2c(mu,nu,jsub,isub)=cmplx(dotprod(1),dotprod(2),kind=dpc)
               end do
             end do
           end do
         end if
       end if

     end do ! End of loop on bands

     if(rfphon==1) then
       do iatom=1,natom
         ia=iatom;if (dyfr_nondiag==0) ia=1
         dyfrnl(1,1,1,iatom,ia)=dyfrnl(1,1,1,iatom,ia)+wtk_k*dyfrnlk(1,iatom)
         dyfrnl(1,2,2,iatom,ia)=dyfrnl(1,2,2,iatom,ia)+wtk_k*dyfrnlk(2,iatom)
         dyfrnl(1,3,3,iatom,ia)=dyfrnl(1,3,3,iatom,ia)+wtk_k*dyfrnlk(3,iatom)
         dyfrnl(1,2,3,iatom,ia)=dyfrnl(1,2,3,iatom,ia)+wtk_k*dyfrnlk(4,iatom)
         dyfrnl(1,1,3,iatom,ia)=dyfrnl(1,1,3,iatom,ia)+wtk_k*dyfrnlk(5,iatom)
         dyfrnl(1,1,2,iatom,ia)=dyfrnl(1,1,2,iatom,ia)+wtk_k*dyfrnlk(6,iatom)
       end do
     end if ! end if rfphon
     if(rfstrs/=0) then
       eltfrnl(:,:)=eltfrnl(:,:)+dtset%wtk(ikpt)*eltfrnlk(:,:)
     end if
     if(need_becfr) then
       becfrnl(:,:,:)=becfrnl(:,:,:)+dtset%wtk(ikpt)*becfrnlk(:,:,:)
     end if
     if(need_piezofr) then
       piezofrnl(:,:)=piezofrnl(:,:)+dtset%wtk(ikpt)*piezofrnlk(:,:)
     end if
!    Increment indexes
     bdtot_index=bdtot_index+nband_k
     if (dtset%mkmem/=0) then
       ibg=ibg+nband_k*dtset%nspinor
       icg=icg+npw_k*dtset%nspinor*nband_k
       ikg=ikg+npw_k
     end if

     ABI_DEALLOCATE(ffnl)
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(ph3d)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylmgr_k)
     ABI_DEALLOCATE(cwavef)
     ABI_DEALLOCATE(kg_k)
     if (need_becfr.or.need_piezofr) then
       ABI_DEALLOCATE(svectout)
     end if
     if (need_piezofr) then
       ABI_DEALLOCATE(ffnl_str)
     end if
     if (need_efmas) then
       ABI_DEALLOCATE(ddkinpw)
       ABI_DEALLOCATE(cg_left)
       ABI_DEALLOCATE(gh2c)
       ABI_DEALLOCATE(gs2c)
     end if

   end do ! End loops on isppol and ikpt
 end do
 if(rfphon==1) then
   ABI_DEALLOCATE(dyfrnlk)
   ABI_DEALLOCATE(enlout_phon)
 end if
 if(rfstrs/=0) then
   ABI_DEALLOCATE(eltfrnlk)
   ABI_DEALLOCATE(enlout_strs)
 end if
 if (need_becfr)  then
   ABI_DEALLOCATE(becfrnlk)
   ABI_DEALLOCATE(enlout_bec1)
   ABI_DEALLOCATE(enlout_bec2)
 end if
 if(need_piezofr)then
   ABI_DEALLOCATE(enlout_piez1)
   ABI_DEALLOCATE(enlout_piez2)
   ABI_DEALLOCATE(piezofrnlk)
 end if
 if(need_efmas) then
   ABI_DEALLOCATE(enlout_efmas)
 end if
 if (psps%usepaw==1) then
   if (need_becfr.or.need_piezofr)  then
     ABI_DEALLOCATE(becij)
   end if
   call pawcprj_free(cwaveprj)
 end if
 ABI_DATATYPE_DEALLOCATE(cwaveprj)

!Fill in lower triangle of matrixes
 if (rfphon==1) then
   do iatom=1,natom
     ia=iatom;if (dyfr_nondiag==0) ia=1
     dyfrnl(1,3,2,iatom,ia)=dyfrnl(1,2,3,iatom,ia)
     dyfrnl(1,3,1,iatom,ia)=dyfrnl(1,1,3,iatom,ia)
     dyfrnl(1,2,1,iatom,ia)=dyfrnl(1,1,2,iatom,ia)
   end do
 end if
 if(rfstrs/=0)then
   do jj=2,6
     do ii=1,jj-1
       eltfrnl(jj,ii)=eltfrnl(ii,jj)
     end do
   end do
 end if

!Parallel case: accumulate (n,k) contributions
 if (xmpi_paral==1) then
   call timab(48,1,tsec)
!  Accumulate dyfrnl
   if(rfphon==1)then
     call xmpi_sum(dyfrnl,spaceworld,ierr)
   end if
!  Accumulate eltfrnl.
   if(rfstrs/=0)then
     call xmpi_sum(eltfrnl,spaceworld,ierr)
   end if
!  Accumulate becfrnl
   if (need_becfr) then
     call xmpi_sum(becfrnl,spaceworld,ierr)
   end if
!  Accumulate piezofrnl
   if (need_piezofr) then
     call xmpi_sum(piezofrnl,spaceworld,ierr)
   end if

!  PAW: accumulate gradients of rhoij
   if (psps%usepaw==1) then
     ABI_ALLOCATE(dimlmn,(natom))
     dimlmn(1:natom)=pawrhoij_tot(1:natom)%cplex_rhoij*pawrhoij_tot(1:natom)%qphase*pawrhoij_tot(1:natom)%lmn2_size
     bufdim=ncpgr*sum(dimlmn)*nsploop
     ABI_ALLOCATE(mpibuf,(bufdim))
     ii=0;mpibuf=zero
     do iatom=1,natom
       do isppol=1,nsploop
         do mu=1,ncpgr
           mpibuf(ii+1:ii+dimlmn(iatom))=pawrhoij_tot(iatom)%grhoij(mu,1:dimlmn(iatom),isppol)
           ii=ii+dimlmn(iatom)
         end do
       end do
     end do
     call xmpi_sum(mpibuf,spaceworld,ierr)
     ii=0
     do iatom=1,natom
       do isppol=1,nsploop
         do mu=1,ncpgr
           pawrhoij_tot(iatom)%grhoij(mu,1:dimlmn(iatom),isppol)=mpibuf(ii+1:ii+dimlmn(iatom))
           ii=ii+dimlmn(iatom)
         end do
       end do
     end do
     ABI_DEALLOCATE(mpibuf)
     ABI_DEALLOCATE(dimlmn)
   end if
   call timab(48,2,tsec)
 end if

!====== PAW: Additional steps
 if (psps%usepaw==1) then

!  Symmetrize rhoij gradients and transfer to cartesian (reciprocal space) coord.
!  This symetrization is necessary in the antiferromagnetic case...
   if (rfphon==1.and.rfstrs==0) then
     option_rhoij=2;option=0
     call pawrhoij_symrhoij(pawrhoij_tot,pawrhoij_tot,option_rhoij,gprimd,indsym,0,natom,dtset%nsym,&
&     psps%ntypat,option,pawang,dtset%pawprtvol,pawtab,rprimd,dtset%symafm,symrec,dtset%typat,&
&     comm_atom=my_comm_atom, mpi_atmtab=my_atmtab)
   else if (rfphon==1.and.rfstrs==1) then
     option_rhoij=23;option=0
     call pawrhoij_symrhoij(pawrhoij_tot,pawrhoij_tot,option_rhoij,gprimd,indsym,0,natom,dtset%nsym,&
&     psps%ntypat,option,pawang,dtset%pawprtvol,pawtab,rprimd,dtset%symafm,symrec,dtset%typat,&
&     comm_atom=my_comm_atom, mpi_atmtab=my_atmtab)
   end if

!  Translate coordinates
   ABI_CHECK(nsploop/=4,'d2frnl: should we mix mx/my/mz when translating coordinates?')
   do iatom=1,natom
     cplx=pawrhoij_tot(iatom)%cplex_rhoij
     do iq=1,pawrhoij_tot(iatom)%qphase
       iq0=0;if (iq==2) iq0=cplx*pawrhoij_tot(iatom)%lmn2_size
       do isppol=1,nsploop
         do klmn=1,pawrhoij_tot(iatom)%lmn2_size
           do ii=1,cplx
             if(rfphon==1.or.rfstrs/=0)then
               grhoij(1:3)=pawrhoij_tot(iatom)%grhoij(shift_rhoij+1:shift_rhoij+3,iq0+cplx*(klmn-1)+ii,isppol)
               do mu=1,3
                 pawrhoij_tot(iatom)%grhoij(shift_rhoij+mu,iq0+cplx*(klmn-1)+ii,isppol)=gprimd(mu,1)*grhoij(1)&
&                  +gprimd(mu,2)*grhoij(2)+gprimd(mu,3)*grhoij(3)
               end do
             end if
             if(rfstrs/=0)then
               call strconv(pawrhoij_tot(iatom)%grhoij(1:6,iq0+cplx*(klmn-1)+ii,isppol),gprimd,&
&                           pawrhoij_tot(iatom)%grhoij(1:6,iq0+cplx*(klmn-1)+ii,isppol))
             end if
           end do
         end do
       end do
     end do
   end do

!  In case of elastic tensor computation, add diagonal contribution:
!     -delta_{alphabeta} rhoi_{ij} to drhoij/d_eps
   if(rfstrs/=0)then
     do iatom=1,natom
       cplx=pawrhoij_tot(iatom)%cplex_rhoij
       do iq=1,pawrhoij_tot(iatom)%qphase
         iq0=0;if (iq==2) iq0=cplx*pawrhoij_tot(iatom)%lmn2_size
         do isppol=1,nsploop
           do nu=1,pawrhoij_tot(iatom)%nrhoijsel
             klmn=pawrhoij_tot(iatom)%rhoijselect(nu)
             do ii=1,cplx
               pawrhoij_tot(iatom)%grhoij(1:3,iq0+cplx*(klmn-1)+ii,isppol)= &
&               pawrhoij_tot(iatom)%grhoij(1:3,iq0+cplx*(klmn-1)+ii,isppol)&
&               -pawrhoij_tot(iatom)%rhoijp(iq0+cplx*(nu-1)+ii,isppol)
             end do
           end do
         end do
       end do
     end do
   end if

!  Add gradients due to Dij derivatives to dynamical matrix/stress tensor
   dimnhat=0;optgr=0;optgr2=0;optstr=0;optstr2=0
   if (rfphon==1) optgr2=1
   if (rfstrs/=0) optstr2=1
   ABI_ALLOCATE(nhat_dum,(1,0))
   call pawgrnl(gs_ham%atindx1,dimnhat,dyfrnl,dyfr_cplex,eltfrnl,dummy,gsqcut,mgfftf,my_natom,natom,&
&   gs_ham%nattyp,nfftf,ngfftf,nhat_dum,dummy,dtset%nspden,dtset%nsym,psps%ntypat,optgr,optgr2,optstr,optstr2,&
&   pawang,pawfgrtab,pawrhoij_tot,pawtab,ph1df,psps,dtset%qptn,rprimd,symrec,dtset%typat,ucvol,vtrial,vxc,xred,&
&   mpi_atmtab=my_atmtab,comm_atom=my_comm_atom)
   ABI_DEALLOCATE(nhat_dum)
 end if !PAW

!The indexing array atindx is used to reestablish the correct order of atoms
 if (rfstrs/=0)then
   ABI_ALLOCATE(elt_work,(6+3*natom,6))
   elt_work(1:6,1:6)=eltfrnl(1:6,1:6)
   do ia=1,natom
     ielt=7+3*(ia-1)
     ieltx=7+3*(gs_ham%atindx(ia)-1)
     elt_work(ielt:ielt+2,1:6)=eltfrnl(ieltx:ieltx+2,1:6)
   end do
   eltfrnl(:,:)=elt_work(:,:)
   ABI_DEALLOCATE(elt_work)
 end if

!Born Effective Charges and PAW:
!1-Re-order atoms -- 2-Add diagonal contribution from rhoij
!3-Multiply by -1 because that the effective charges
 !  are minus the second derivatives of the energy
 if (need_becfr) then
   ABI_ALLOCATE(becfrnl_tmp,(3,natom,3))
   becfrnl_tmp=-becfrnl
   do ia=1,natom         ! Atom (sorted by type)
     iatom=gs_ham%atindx1(ia)   ! Atom (not sorted)
     itypat=dtset%typat(iatom)
     do ii=1,3           ! Direction of electric field
       do jj=1,3         ! Direction of atom
         becfrnl(jj,iatom,ii)=becfrnl_tmp(jj,ia,ii)
       end do
     end do
   end do
   ABI_DEALLOCATE(becfrnl_tmp)
 end if

!Piezoelectric Tensor
!-Multiply by -1 because that the piezoelectric tensor
!  are minus the second derivatives of the energy
 if (need_piezofr) then
   piezofrnl=-piezofrnl
 end if

!Close the ddk files
 do ii=1,3
   call ddkfiles(ii)%close()
 end do

!Release now useless memory
 if (psps%usepaw==1) then
   do iatom=1,natom
     ABI_DEALLOCATE(pawrhoij_tot(iatom)%grhoij)
     pawrhoij_tot(iatom)%ngrhoij=0
   end do
   if (paral_atom) then
     call pawrhoij_free(pawrhoij_tot)
     ABI_DATATYPE_DEALLOCATE(pawrhoij_tot)
   end if
 end if
 call gs_ham%free()
 call timab(159,2,tsec)

 write(msg,'(3a)')ch10,' ==> Calculation of the frozen part of the second order derivative done',ch10
 call wrtout(std_out,msg,'COLL')

 DBG_EXIT("COLL")

end subroutine d2frnl
!!***

end module m_d2frnl
!!***
