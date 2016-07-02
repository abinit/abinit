!{\src2tex{textfont=tt}}
!!****f* ABINIT/posdoppler
!! NAME
!! posdoppler
!!
!! FUNCTION
!! Calculate the momentum distribution annihilating electrons-positron (Doppler broadening)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (JW,GJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!!      set_mpi_enreg_fft,simp_gen,sphereboundary,symrhoij,unset_mpi_enreg_fft
!!      wffclose,wffopen,wrtout,xderivewrecend,xderivewrecinit,xderivewrite
!!      xmoveoff,xmpi_bcast,xmpi_recv,xmpi_send,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

!Macro to go from row-column indexing to combined indexing
#define RCC(glmn,hlmn) max(glmn,hlmn)*(max(glmn,hlmn)-1)/2+min(glmn,hlmn)
!Macro to go from l,m angular momentum indexing to combined indexing
#define LMC(lval,mval) lval*lval+lval+mval+1

subroutine posdoppler(cg,cprj,Crystal,dimcprj,dtfil,dtset,electronpositron,&
&                     filpsp,kg,mcg,mcprj,mpi_enreg,my_natom,&
&                     n3xccc,nfft,ngfft,nhat,npwarr,occ,pawang,pawrad,&
&                     pawrhoij,pawtab,rhor,xccc3d)


 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_bandfft_kpt
 use m_wffile
 use m_electronpositron

 use m_pawcprj, only : pawcprj_type, pawcprj_alloc, pawcprj_get, pawcprj_mpi_send, &
&                      pawcprj_mpi_recv, pawcprj_free, pawcprj_copy, pawcprj_bcast
 use m_pawang,  only : pawang_type, realgaunt
 use m_pawrad,  only : pawrad_type, simp_gen
 use m_pawtab,  only : pawtab_type
 use m_pawrhoij,only : pawrhoij_type, pawrhoij_alloc, pawrhoij_free,&
&                      pawrhoij_nullify, pawrhoij_gather, symrhoij
 use m_pawxc,   only : pawxcsum
 use m_paw_sphharm, only : initylmr
 use m_pawpsp,  only : pawpsp_read_corewf
 use m_crystal, only : crystal_t
 use m_mpinfo,  only : ptabs_fourdp,set_mpi_enreg_fft,unset_mpi_enreg_fft,destroy_mpi_enreg
 use m_io_tools,only : open_file,close_unit,get_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'posdoppler'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_xc_lowlevel
 use interfaces_51_manage_mpi
 use interfaces_52_fft_mpi_noabirule
 use interfaces_53_ffts
 use interfaces_56_xc
 use interfaces_65_paw
 use interfaces_66_wfs
 use interfaces_67_common, except_this_one => posdoppler
!End of the abilint section

 implicit none

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
             pbn(:) = pcart(:)/pnorm ! unit vector

             if (pnorm < tol12) then
               pbn(:) = zero
               ylmp(:) = zero
               ylmp(1) = 1.d0/sqrt(four_pi)
             else
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
     nspden_rhoij=dtset%nspden;if (dtset%pawspnorb>0.and.dtset%nspinor==2) nspden_rhoij=4
     call pawrhoij_alloc(pawrhoij_dop_el,dtset%pawcpxocc,nspden_rhoij,&
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
&         n4,n5,n6,option,mpi_enreg%paral_kgb,tim_fourwf,weight_pos,weight_pos,&
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
&                       n4,n5,n6,option,mpi_enreg%paral_kgb,tim_fourwf,weight,weight,&
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
!                            If not, call symrhoij with nsym=1
                             call symrhoij(pawrhoij_dop_el,pawrhoij_dop_el,1,Crystal%gprimd,Crystal%indsym,0,dtset%natom,&
&                             Crystal%nsym,dtset%ntypat,1,pawang,-10001,pawtab,Crystal%rprimd,Crystal%symafm,&
&                             Crystal%symrec,dtset%typat)
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
                         call fourdp(cplex,rho_contrib_g,rho_contrib,-1,mpi_enreg,nfft,ngfft,&
&                         mpi_enreg%paral_kgb,tim_fourdp)

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
