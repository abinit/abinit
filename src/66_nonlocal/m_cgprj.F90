!!****m* ABINIT/m_cgprj
!! NAME
!!  m_cgprj
!!
!! FUNCTION
!!   Routines to compute <Proj_i|Cnk> with |Cnk> expressed in reciprocal space.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (MT)
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

module m_cgprj

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi

 use defs_abitypes, only : MPI_type
 use defs_datatypes, only : pseudopotential_type
 use m_kg,       only : ph1d3d, mkkpg
 use m_geometry, only : strconv
 use m_mkffnl,   only : mkffnl
 use m_mpinfo,   only : proc_distrb_cycle
 use m_pawcprj,  only : pawcprj_type, pawcprj_alloc, pawcprj_put, pawcprj_free, &
                        pawcprj_set_zero, pawcprj_mpi_sum
 use m_opernla_ylm, only : opernla_ylm

 implicit none

 private
!!***

 public :: getcprj
 public :: ctocprj
!!***

contains
!!***

!!****f* ABINIT/getcprj
!! NAME
!! getcprj
!!
!! FUNCTION
!!  Compute <Proj_i|Cnk> for one wave function |Cnk> expressed in reciprocal space.
!!  Compute also derivatives of <Proj_i|Cnk>.
!!  |Proj_i> are non-local projectors (for each atom and each l,m,n)
!!
!! INPUTS
!!  choice=chooses possible output:
!!    In addition to projected wave function:
!!    choice=1 => nothing else
!!          =2 => 1st gradients with respect to atomic position(s)
!!          =3 => 1st gradients with respect to strain(s)
!!          =23=> 1st gradients with respect to strain(s) and atm pos
!!          =4 => 2nd derivatives with respect to atomic pos.
!!          =24=> 1st and 2nd derivatives with respect to atomic pos.
!!          =5 => 1st gradients with respect to k wavevector
!!          =6 => 2nd derivatives with respect to strain and atm. pos.
!!  cpopt=1 if <Proj_i|Cnk> are already in memory; see below (side effects).
!!  cwavef(2,nspinor*npw_k)=input cmplx wavefunction coefficients <G|Cnk>
!!  ffnl(npw_k,dimffnl,lmnmax,ntypat)=nonlocal form factors to be used for the application of the nl operator
!!  idir=direction of the derivative, i.e. dir. of - atom to be moved  in the case choice=2
!!                                                 - strain component  in the case choice=3
!!                                                 - k point direction in the case choice=5
!!       Compatible only with choice=2,3,5; if idir=0, all derivatives are computed
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=lmn
!!  istwf_k=option parameter that describes the storage of wfs
!!  kg_k(3,npw_k)=reduced planewave coordinates
!!  kpg(npw_k,npk)=(k+G) components and related data
!!  kpoint(3)=k point in terms of recip. translations
!!  lmnmax=max. number of (l,m,n) components over all types of atoms
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  nattyp(ntypat)=number of atoms of each type
!!  ngfft(18)=contain all needed information about 3D FFT, see ~ABINIT/Infos/vargs.htm#ngfft
!!  nloalg(3)=governs the choice of the algorithm for nonlocal operator
!!  npw_k=number of planewaves for given k point
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  ntypat=number of types of atoms in unit cell
!!  phkxred(2,natom)=phase factors exp(2 pi kpoint.xred)
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information
!!  ph3d(2,npw_k,natom)=3D structure factors, for each atom and plane wave
!!                      only used if nloalg(2)>0
!!  ucvol= unit cell volume
!!  useylm=governs the way the nonlocal operator is to be applied
!!
!! SIDE EFFECTS
!!  cwaveprj(natom,nspinor) <type(pawcprj_type)>=projected input wave function <Proj_i|Cnk> with all NL projectors
!!                                (and derivatives)
!!                                if cpopt=1 the projected scalars have been already been computed and
!!                                           only derivatives are computed here
!!                                if cpopt=0 the projected scalars and derivatives are computed here
!!
!! TODO
!!  Spin-orbit
!!
!! PARENTS
!!      cgwf,ctocprj,debug_tools,dfpt_accrho,dfpt_nstpaw,ks_ddiago,m_wfd
!!      rf2_init
!!
!! CHILDREN
!!      mkkpg,opernla_ylm,ph1d3d
!!
!! SOURCE

 subroutine getcprj(choice,cpopt,cwavef,cwaveprj,ffnl,&
&                   idir,indlmn,istwf_k,kg_k,kpg,kpoint,lmnmax,mgfft,mpi_enreg,&
&                   natom,nattyp,ngfft,nloalg,npw_k,nspinor,ntypat,&
&                   phkxred,ph1d,ph3d,ucvol,useylm)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: choice,cpopt,idir,istwf_k,lmnmax
 integer,intent(in) :: mgfft,natom,npw_k,nspinor,ntypat,useylm
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),kg_k(3,npw_k),nattyp(ntypat)
 integer,intent(in) :: ngfft(18),nloalg(3)
 real(dp),intent(in) :: cwavef(2,npw_k*nspinor)
 real(dp),intent(in),target :: ffnl(:,:,:,:),kpg(:,:),ph3d(:,:,:)
 real(dp),intent(in) :: kpoint(3),ph1d(2,3*(2*mgfft+1)*natom),phkxred(2,natom)
 type(pawcprj_type),intent(inout) :: cwaveprj(natom,nspinor)

!Local variables-------------------------------
!scalars
 integer :: choice_,cplex,dimffnl,ia,ia1,ia2,ia3,ia4,iatm,ic,ii,ilmn,ishift,ispinor,itypat
 integer :: jc,matblk,mincat,nd2gxdt,ndgxdt,nincat,nkpg,nkpg_,nlmn,signs
!arrays
 integer,allocatable :: cplex_dgxdt(:),cplex_d2gxdt(:),indlmn_typ(:,:)
 real(dp),allocatable :: d2gxdt(:,:,:,:,:),dgxdt(:,:,:,:,:),ffnl_typ(:,:,:)
 real(dp),allocatable :: gx(:,:,:,:)
 real(dp), pointer :: kpg_(:,:),ph3d_(:,:,:)

! *********************************************************************

 DBG_ENTER('COLL')

!Nothing to do in that case
 if (cpopt==1.and.choice==1) return

!Not available for useylm=0
 if (useylm==0) then
   MSG_ERROR('Not available for useylm=0 !')
 end if

!Error on bad choice
 if ((choice<1.or.choice>6).and.choice/=23.and.choice/=24) then
   MSG_BUG('Does not presently support this choice !')
 end if

!Error on bad idir
 if (idir>0.and.choice/=2.and.choice/=3.and.choice/=5) then
   MSG_BUG('Does not support idir>0 for this choice')
 end if

!Error on sizes
 nkpg=size(kpg,2)
 if (nkpg>0) then
   if( (choice==2.and.nkpg<3) .or. &
&   ((choice==4.or.choice==24).and.nkpg<9) .or. &
&   ((choice==6.or.choice==3.or.choice==23).and.nkpg<3) ) then
     MSG_BUG('Incorrect size for kpg array !')
   end if
 end if
 if (size(ffnl,1)/=npw_k.or.size(ffnl,3)/=lmnmax) then
   MSG_BUG('Incorrect size for ffnl!')
 end if
 if (size(ph3d)>0) then
   if (size(ph3d,2)/=npw_k) then
     MSG_BUG('Incorrect size for ph3d!')
   end if
 end if

!Define dimensions of projected scalars
 dimffnl=size(ffnl,2)
 ndgxdt=0;nd2gxdt=0
 if (idir==0) then
   if (choice==2) ndgxdt=3
   if (choice==3) ndgxdt=6
   if (choice==23) ndgxdt=9
   if (choice==4) nd2gxdt=6
   if (choice==24) then
     ndgxdt=3;nd2gxdt=6
   end if
   if (choice==5) ndgxdt=3
   if (choice==6) then
     ndgxdt=9;nd2gxdt=54
   end if
 else
   ndgxdt=1
 end if
 if(cwaveprj(1,1)%ncpgr<ndgxdt+nd2gxdt) then
   MSG_BUG('Incorrect size for ncpgr')
 end if
!Eventually re-compute (k+G) vectors (and related data)
 if (nkpg==0) then
   nkpg_=0
   if (choice==4.or.choice==24) nkpg_=9
   if (choice==2.or.choice==3.or.choice==23) nkpg_=3
   ABI_ALLOCATE(kpg_,(npw_k,nkpg_))
   if (nkpg_>0) then
     call mkkpg(kg_k,kpg_,kpoint,nkpg_,npw_k)
   end if
 else
   nkpg_=nkpg
   kpg_ => kpg
 end if

!Some other dims
 mincat=min(NLO_MINCAT,maxval(nattyp))
 cplex=2;if (istwf_k>1) cplex=1
 choice_=choice;if (cpopt==1) choice_=-choice
 signs=1;if (idir>0) signs=2
!Eventually allocate temporary array for ph3d
 if (nloalg(2)<=0) then
   matblk=mincat
   ABI_ALLOCATE(ph3d_,(2,npw_k,matblk))
 else
   matblk=size(ph3d,3)
   ph3d_ => ph3d
 end if

!Loop over atom types
 ia1=1;iatm=0
 do itypat=1,ntypat
   ia2=ia1+nattyp(itypat)-1;if (ia2<ia1) cycle
   nlmn=count(indlmn(3,:,itypat)>0)

!  Retrieve some data for this type of atom
   ABI_ALLOCATE(indlmn_typ,(6,nlmn))
   ABI_ALLOCATE(ffnl_typ,(npw_k,dimffnl,nlmn))
   indlmn_typ(:,1:nlmn)=indlmn(:,1:nlmn,itypat)
   ffnl_typ(:,:,1:nlmn)=ffnl(:,:,1:nlmn,itypat)

!  Loop on blocks of atoms inside type
   do ia3=ia1,ia2,mincat
     ia4=min(ia2,ia3+mincat-1);nincat=ia4-ia3+1
!     Prepare the phase factors if they were not already computed
     if (nloalg(2)<=0) then
       call ph1d3d(ia3,ia4,kg_k,matblk,natom,npw_k,ngfft(1),ngfft(2),ngfft(3),&
&       phkxred,ph1d,ph3d_)
     end if

!    Allocate memory for projected scalars
     ABI_ALLOCATE(gx,(cplex,nlmn,nincat,nspinor))
     ABI_ALLOCATE(dgxdt,(cplex,ndgxdt,nlmn,nincat,nspinor))
     ABI_ALLOCATE(d2gxdt,(cplex,nd2gxdt,nlmn,nincat,nspinor))
     ABI_ALLOCATE(cplex_dgxdt,(ndgxdt))
     ABI_ALLOCATE(cplex_d2gxdt,(nd2gxdt))

!    Retrieve eventually <p_i|c> coeffs
     if (cpopt==1) then
       do ispinor=1,nspinor
         do ia=1,nincat
           gx(1:cplex,1:nlmn,ia,ispinor)=cwaveprj(iatm+ia,ispinor)%cp(1:cplex,1:nlmn)
         end do
       end do
     end if

!    Compute <p_i|c> scalars (and derivatives) for this block of atoms
     call opernla_ylm(choice_,cplex,cplex_dgxdt,cplex_d2gxdt,dimffnl,d2gxdt,dgxdt,ffnl_typ,gx,&
&     ia3,idir,indlmn_typ,istwf_k,kpg_,matblk,mpi_enreg,nd2gxdt,ndgxdt,nincat,nkpg_,nlmn,&
&     nloalg,npw_k,nspinor,ph3d_,signs,ucvol,cwavef)

!    Transfer result to output variable cwaveprj
     if (cpopt==0) then
       do ispinor=1,nspinor
         do ia=1,nincat
           cwaveprj(iatm+ia,ispinor)%nlmn=nlmn
           cwaveprj(iatm+ia,ispinor)%cp(1:cplex,1:nlmn)=gx(1:cplex,1:nlmn,ia,ispinor)
           if(cplex==1) cwaveprj(iatm+ia,ispinor)%cp(2,1:nlmn)=zero
         end do
       end do
     end if
     if (cpopt>=0.and.choice>1) then
       ishift=0
       if ((idir>0).and.(cwaveprj(1,1)%ncpgr>ndgxdt)) ishift=idir-1
       if(cplex==2)then
         do ispinor=1,nspinor
           do ia=1,nincat
!             cwaveprj(iatm+ia,ispinor)%ncpgr=ndgxdt+nd2gxdt
             if (ndgxdt>0) cwaveprj(iatm+ia,ispinor)%dcp(1:2,1+ishift:ndgxdt+ishift,1:nlmn)=&
&             dgxdt(1:2,1:ndgxdt,1:nlmn,ia,ispinor)
             if (nd2gxdt>0)cwaveprj(iatm+ia,ispinor)%dcp(1:2,ndgxdt+1+ishift:ndgxdt+nd2gxdt+ishift,1:nlmn)=&
&             d2gxdt(1:2,1:nd2gxdt,1:nlmn,ia,ispinor)
           end do
         end do
       else
!        cplex_dgxdt(i)  = 1 if dgxdt(1,i,:,:)  is real, 2 if it is pure imaginary
!        cplex_d2gxdt(i) = 1 if d2gxdt(1,i,:,:) is real, 2 if it is pure imaginary
         do ispinor=1,nspinor
           do ia=1,nincat
!             cwaveprj(iatm+ia,ispinor)%ncpgr=ndgxdt+nd2gxdt
             if (ndgxdt>0) then
               do ilmn =1,nlmn
                 do ii = 1,ndgxdt
                   ic = cplex_dgxdt(ii) ; jc = 3 - ic
                   cwaveprj(iatm+ia,ispinor)%dcp(ic,ii+ishift,ilmn)=dgxdt(1,ii,ilmn,ia,ispinor)
                   cwaveprj(iatm+ia,ispinor)%dcp(jc,ii+ishift,ilmn)=zero
                 end do
               end do
             end if
             if (nd2gxdt>0) then
               do ilmn =1,nlmn
                 do ii = 1,nd2gxdt
                   ic = cplex_d2gxdt(ii) ; jc = 3 - ic
                   cwaveprj(iatm+ia,ispinor)%dcp(ic,ndgxdt+ii+ishift,ilmn)=d2gxdt(1,ii,ilmn,ia,ispinor)
                   cwaveprj(iatm+ia,ispinor)%dcp(jc,ndgxdt+ii+ishift,ilmn)=zero
                 end do
               end do
             end if
           end do
         end do
       end if
     end if

!    End loop inside block of atoms
     iatm=iatm+nincat
     ABI_DEALLOCATE(gx)
     ABI_DEALLOCATE(dgxdt)
     ABI_DEALLOCATE(d2gxdt)
     ABI_DEALLOCATE(cplex_dgxdt)
     ABI_DEALLOCATE(cplex_d2gxdt)
   end do

!  End loop over atom types
   ia1=ia2+1
   ABI_DEALLOCATE(indlmn_typ)
   ABI_DEALLOCATE(ffnl_typ)
 end do

 if (nkpg==0) then
   ABI_DEALLOCATE(kpg_)
 end if
 if (nloalg(2)<=0) then
   ABI_DEALLOCATE(ph3d_)
 end if

 DBG_EXIT('COLL')

 end subroutine getcprj
!!***

!!****f* ABINIT/ctocprj
!! NAME
!! ctocprj
!!
!! FUNCTION
!!  Compute all <Proj_i|Cnk> for every wave function |Cnk> expressed in reciprocal space.
!!  |Proj_i> are non-local projectors (for each atom and each l,m,n)
!!  Can also compute derivatives of <Proj_i|Cnk> wrt to several parameters
!!
!! INPUTS
!!  atindx(natom)=index table for atoms
!!  cg(2,mcg)=planewave coefficients of wavefunctions
!!  choice: chooses derivatives to compute:
!!          =1 => no derivatives
!!          =2 => 1st derivatives with respect to atomic position(s)
!!          =3 => 1st derivatives with respect to strain(s)
!!          =23=> 1st derivatives with respect to strain(s) and atm pos
!!          =4 => 2nd derivatives with respect to atomic pos.
!!          =24=> 1st and 2nd derivatives with respect to atomic pos.
!!          =5 => derivatives with respect to k wavevector
!!          =6 => 2nd derivatives with respect to strain and atm. pos.
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  iatom= if <=0, cprj=<p_i|Cnk> are computed for all atoms 1...natom
!!         if  >0  cprj=<p_i|Cnk> are computed only for atom with index iatom
!!  idir=direction of the derivative, i.e. dir. of - atom to be moved  in the case choice=2
!!                                                 - strain component  in the case choice=3
!!                                                 - k point direction in the case choice=5
!!       Compatible only with choice=2,3,5; if idir=0, all derivatives are computed
!!  iorder_cprj=0 if output cprj=<p_i|Cnk> are sorted by atom type
!!                (first all elements of atom type 1, followed by those of atom type 2 and so on).
!!              1 if output cprj=<p_i|Cnk> are sorted according to
!!                the variable typat in the main input file
!!  istwfk(nkpt)=option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates
!!  kpt(3,nkpt)=reduced coordinates of k points.
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband_mem*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mgfft=maximum size of 1D FFTs
!!  mkmem=number of k points treated by this node.
!!  mpi_enreg=information about MPI parallelization
!!  mpsang=1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw
!!  natom=number of atoms in cell
!!  nattyp(ntypat)= # atoms of each type
!!  nband(nkpt*nsppol)=number of bands at this k point for that spin polarization
!TODO : distribute cprj over bands as well
!!  mband_mem=max number of bands for this processor (in case of band parallelism)
!!  ncprj=1st dim. of cprj array (natom if iatom<=0, 1 if iatom>0)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~ABINIT/Infos/vargs.htm#ngfft
!!  nkpt=number of k points
!!  nloalg(3)=governs the choice of the algorithm for nonlocal operator
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell
!!  paral_kgb= 1 if kpt-band-FFT is activated
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rmet(3,3)=real space metric (bohr**2)
!!  tim_ctocprj=timing code of the calling routine
!!  typat(natom)= types of atoms
!!  uncp=unit number for <P_lmn|Cnk> data (if used)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang)=real spherical harmonics for each G and k point
!!  ylmgr(npw*mkmem,nylmgr,mpsang*mpsang*useylmgr)=gradients of real spherical harmonics wrt (k+G)
!!!
!! OUTPUT
!!  cprj(ncprj,mcprj) <type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk> with NL projectors
!!                                       Usually ncprj=natom
!!
!! PARENTS
!!      dfpt_looppert,extrapwf,forstr,scfcv,update_e_field_vars,vtorho
!!
!! CHILDREN
!!      getcprj,mkffnl,mkkpg,pawcprj_alloc,pawcprj_free,pawcprj_mpi_sum
!!      pawcprj_put,pawcprj_set_zero,ph1d3d,strconv,xmpi_allgather
!!      xmpi_allgatherv,xmpi_alltoallv
!!
!! SOURCE

 subroutine ctocprj(atindx,cg,choice,cprj,gmet,gprimd,iatom,idir,&
& iorder_cprj,istwfk,kg,kpt,mcg,mcprj,mgfft,mkmem,mpi_enreg,mpsang,&
& mpw,natom,nattyp,nband,ncprj,ngfft,nkpt,nloalg,npwarr,nspinor,&
& nsppol,ntypat,paral_kgb,ph1d,psps,rmet,typat,ucvol,uncp,xred,ylm,ylmgr)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: choice,iatom,idir,iorder_cprj,mcg,mcprj,mgfft,mkmem,mpsang,mpw
 integer,intent(in) :: natom,ncprj,nkpt,nspinor,nsppol,ntypat,paral_kgb,uncp
!TODO : distribute cprj over bands as well
! integer,intent(in) :: mband_mem
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
 type(pseudopotential_type),target,intent(in) :: psps
!arrays
 integer,intent(in) :: istwfk(nkpt),nband(nkpt*nsppol)
 integer,intent(in) :: ngfft(18),nloalg(3),npwarr(nkpt),kg(3,mpw*mkmem),typat(natom)
 integer,intent(in),target :: atindx(natom),nattyp(ntypat)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),kpt(3,nkpt),rmet(3,3)
 real(dp),intent(in) :: xred(3,natom),ylm(:,:),ylmgr(:,:,:)
 real(dp),intent(in),target :: ph1d(2,3*(2*mgfft+1)*natom)
 type(pawcprj_type),intent(inout) :: cprj(ncprj,mcprj)

!Local variables-------------------------------
!scalars
 integer :: blocksz,cg_bandpp,counter,cpopt,cprj_bandpp,dimffnl,ia,iatm,iatom1,iatom2
 integer :: iband_max,iband_min,iband_start,ibg,ibgb,iblockbd,ibp,icg,icgb,icp1,icp2
 integer :: ider,idir0,iend,ierr,ig,ii,ikg,ikpt,ilm,ipw,isize,isppol,istart,istwf_k,itypat,iwf1,iwf2,jdir
 integer :: matblk,mband_cprj,me_distrb,my_nspinor,n1,n1_2p1,n2,n2_2p1,n3,n3_2p1,kk,nlmn
 integer :: nband_k,nband_cprj_k,nblockbd,ncpgr,nkpg,npband_bandfft,npws,npw_k,npw_nk,ntypat0
 integer :: shift1,shift1b,shift2,shift2b,shift3,shift3b
 integer :: spaceComm,spaceComm_band,spaceComm_fft,useylmgr
 logical :: cg_band_distributed,cprj_band_distributed,one_atom
 real(dp) :: arg
 character(len=500) :: msg
!arrays
 integer,allocatable :: bufsize(:),bufsize_wf(:),bufdisp(:),bufdisp_wf(:)
 integer,allocatable :: dimlmn(:),kg_k(:,:),kg_k_loc(:,:)
 integer,allocatable :: npw_block(:),npw_disp(:)
 integer,pointer :: atindx_atm(:),indlmn_atm(:,:,:),nattyp_atm(:),pspso_atm(:)
 real(dp) :: kpoint(3),work(6)
 real(dp),allocatable :: cwavef(:,:),cwavef_tmp(:,:)
 real(dp),allocatable :: ffnl(:,:,:,:),ffnl_npw(:,:,:,:),ffnl_tmp(:,:,:,:),ffnl_tmp_npw(:,:,:,:)
 real(dp),allocatable :: kpg_k(:,:)
 real(dp),allocatable :: ph3d(:,:,:),ph3d_npw(:,:,:),ph3d_tmp(:,:,:),ph3d_tmp_npw(:,:,:)
 real(dp),allocatable :: phkxred(:,:),ylm_k(:,:),ylmgr_k(:,:,:)
 real(dp),ABI_CONTIGUOUS pointer :: ekb_atm(:,:),ffspl_atm(:,:,:,:),ph1d_atm(:,:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)

! *********************************************************************

 DBG_ENTER('COLL')

!Nothing to do if current MPI process does treat kpoints or plane-waves
 if (mcg==0.or.mcprj==0) return

!Preliminary tests
 if (psps%useylm==0) then
   msg='Not available for useylm=0!'
   MSG_ERROR(msg)
 end if
 if ((choice<1.or.choice>6).and.choice/=23.and.choice/=24) then
   msg='Bad choice!'
   MSG_BUG(msg)
 end if
 if (idir>0.and.choice/=2.and.choice/=3.and.choice/=5) then
   msg='Does not support idir>0 for that choice!'
   MSG_BUG(msg)
 end if
 if (size(ylm)/=mpw*mkmem*mpsang*mpsang) then
   msg='Wrong size for Ylm!'
   MSG_BUG(msg)
 end if
 useylmgr=merge(0,1,(size(ylmgr)==0))
 if (useylmgr==0.and.(choice==3.or.choice==5.or.choice==23)) then
   msg=' Ylm gradients have to be in memory for choice=3, 5, or 23!'
   MSG_BUG(msg)
 end if

!Init parallelism
 if (paral_kgb==1) then
   me_distrb=mpi_enreg%me_kpt
   spaceComm=mpi_enreg%comm_kpt
   spaceComm_fft=mpi_enreg%comm_fft
   npband_bandfft=mpi_enreg%nproc_band
   cg_bandpp=mpi_enreg%bandpp
   cprj_bandpp=mpi_enreg%bandpp
   spaceComm_band=mpi_enreg%comm_band
   cg_band_distributed=.true.
   cprj_band_distributed=(mcprj/=mcg/mpw)
 else
   me_distrb=mpi_enreg%me_kpt
   spaceComm=mpi_enreg%comm_cell
   spaceComm_fft=xmpi_comm_self
   npband_bandfft=1
   cg_bandpp=1
   cprj_bandpp=1
   cg_band_distributed=.false.
   cprj_band_distributed=.false.
   spaceComm_band=xmpi_comm_self
   if (mpi_enreg%paralbd==1) then
     cg_band_distributed=.true.
     cprj_band_distributed=.true.
     spaceComm_band=mpi_enreg%comm_band
   end if
 end if
 if (cg_bandpp/=cprj_bandpp) then
   MSG_BUG('cg_bandpp must be equal to cprj_bandpp !')
 end if

!Manage parallelization over spinors
 my_nspinor=max(1,nspinor/mpi_enreg%nproc_spinor)
!Check sizes for cprj (distribution is tricky)
 one_atom=(iatom>0)
 if (one_atom.and.ncprj/=1) then
   MSG_BUG('Bad value for ncprj dimension (should be 1) !')
 end if
 if (.not.one_atom.and.ncprj/=natom) then
   MSG_BUG('Bad value for ncprj dimension (should be natom) !')
 end if

!Initialize some variables
 mband_cprj=mcprj/(my_nspinor*mkmem*nsppol)
 n1=ngfft(1);n2=ngfft(2);n3=ngfft(3)
 n1_2p1=2*n1+1;n2_2p1=2*n2+1;n3_2p1=2*n3+1
 ibg=0;icg=0;cpopt=0
 ider=0;idir0=0;istart=idir;iend=idir
 if (choice==3.or.choice==5.or.choice==23) ider=1
 if (idir>0) then
   if (choice==3) idir0=-idir
   if (choice==5) idir0=idir
 else
!   if (choice==23) idir0=-7
   if (choice==3) idir0=-7
   if (choice==5) idir0=4
 end if
 if (idir0==0.or.idir0==4) then
   dimffnl=1+3*ider
 else if (idir0/=-7) then
   dimffnl=1+ider
 else
   dimffnl=1+6*ider
   if(choice==3)then
     istart=ider;iend=6*ider
   end if
 end if
 nkpg=0
 if (choice==3.or.choice==2.or.choice==23) nkpg=3*nloalg(3)
 if (choice==4.or.choice==24) nkpg=9*nloalg(3)

!Set number of gradients for <p_i|Cnk>
 ncpgr=0
 if (idir==0) then
   if (choice==2) ncpgr=3
   if (choice==3) ncpgr=6
   if (choice==23)ncpgr=9
   if (choice==4) ncpgr=6
   if (choice==24)ncpgr=9
   if (choice==5) ncpgr=3
   if (choice==6) ncpgr=63
 else
   ncpgr=1
 end if
!Test cprj gradients dimension (just to be sure)
 if (cprj(1,1)%ncpgr/=ncpgr) then
   MSG_BUG('cprj are badly allocated !')
 end if


!Extract data for treated atom(s)
 if (one_atom) then
   iatom1=iatom;iatom2=iatom
   ntypat0=1;itypat=typat(iatom)
   ABI_ALLOCATE(nattyp_atm,(ntypat0))
   nattyp_atm(1)=1
   ABI_ALLOCATE(atindx_atm,(ntypat0))
   atindx_atm(1)=atindx(iatom)
   ABI_ALLOCATE(ph1d_atm,(2,(n1_2p1+n2_2p1+n3_2p1)*ntypat0))
   shift1=(atindx(iatom)-1)*n1_2p1
   shift2=(atindx(iatom)-1)*n2_2p1+natom*n1_2p1
   shift3=(atindx(iatom)-1)*n3_2p1+natom*(n1_2p1+n2_2p1)
   shift1b=0;shift2b=n1_2p1;shift3b=n1_2p1+n2_2p1
   ph1d_atm(:,shift1b+1:shift1b+n1_2p1)=ph1d(:,shift1+1:shift1+n1_2p1)
   ph1d_atm(:,shift2b+1:shift2b+n2_2p1)=ph1d(:,shift2+1:shift2+n2_2p1)
   ph1d_atm(:,shift3b+1:shift3b+n3_2p1)=ph1d(:,shift3+1:shift3+n3_2p1)
   ABI_ALLOCATE(ekb_atm,(psps%dimekb,ntypat0))
   ABI_ALLOCATE(indlmn_atm,(6,psps%lmnmax,ntypat0))
   ABI_ALLOCATE(ffspl_atm,(psps%mqgrid_ff,2,psps%lnmax,ntypat0))
   ABI_ALLOCATE(pspso_atm,(ntypat0))
   ekb_atm(:,1)=psps%ekb(:,itypat)
   indlmn_atm(:,:,1)=psps%indlmn(:,:,itypat)
   ffspl_atm(:,:,:,1)=psps%ffspl(:,:,:,itypat)
   pspso_atm(1)=psps%pspso(itypat)
 else
   iatom1=1;iatom2=natom
   ntypat0=ntypat
   atindx_atm => atindx
   nattyp_atm => nattyp
   ph1d_atm => ph1d
   ekb_atm => psps%ekb
   indlmn_atm => psps%indlmn
   ffspl_atm => psps%ffspl
   pspso_atm => psps%pspso
 end if

!Dimensioning and allocation of <p_i|Cnk>
 ABI_ALLOCATE(dimlmn,(ncprj))
 dimlmn=0  ! Type-sorted cprj
 if (one_atom) then
   itypat=typat(iatom)
   dimlmn(ia+1:ia+nattyp(itypat))=count(indlmn_atm(3,:,itypat)>0)
 else
   ia=0
   do itypat=1,ntypat0
     dimlmn(ia+1:ia+nattyp(itypat))=count(indlmn_atm(3,:,itypat)>0)
     ia=ia+nattyp(itypat)
   end do
 end if
 ABI_DATATYPE_ALLOCATE(cwaveprj,(ncprj,my_nspinor*cprj_bandpp))
 call pawcprj_alloc(cwaveprj,ncpgr,dimlmn)

!Additional statements if band-fft parallelism
 if (npband_bandfft>1) then
   ABI_ALLOCATE(npw_block,(npband_bandfft))
   ABI_ALLOCATE(npw_disp,(npband_bandfft))
   ABI_ALLOCATE(bufsize,(npband_bandfft*cg_bandpp))
   ABI_ALLOCATE(bufdisp,(npband_bandfft*cg_bandpp))
   ABI_ALLOCATE(bufsize_wf,(npband_bandfft*cg_bandpp))
   ABI_ALLOCATE(bufdisp_wf,(npband_bandfft*cg_bandpp))
 end if

!Set output datastructure to zero
 call pawcprj_set_zero(cprj)

!LOOP OVER SPINS
 do isppol=1,nsppol
   ikg=0

!  BIG FAT k POINT LOOP
   do ikpt=1,nkpt
     counter=100*ikpt+isppol

!    Select k point to be treated by this proc
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) cycle

!    Retrieve k-point
     kpoint(:)=kpt(:,ikpt)
     istwf_k=istwfk(ikpt)

!    Retrieve number of plane waves
     npw_k=npwarr(ikpt)
     if (npband_bandfft>1) then
!      Special treatment for band-fft //
       call xmpi_allgather(npw_k,npw_block,spaceComm_band,ierr)
       npw_nk=sum(npw_block);npw_disp(1)=0
       do ii=2,npband_bandfft
         npw_disp(ii)=npw_disp(ii-1)+npw_block(ii-1)
       end do
     else
       npw_nk=npw_k
     end if

!    Retrieve (k+G) points and spherical harmonics
     ABI_ALLOCATE(ylm_k,(npw_k,mpsang*mpsang))
     ABI_ALLOCATE(ylmgr_k,(npw_k,3,mpsang*mpsang*useylmgr))
     ABI_ALLOCATE(kg_k,(3,npw_nk))
     if (npband_bandfft>1) then
!      Special treatment for band-fft //
       ABI_ALLOCATE(kg_k_loc,(3,npw_k))
       kg_k_loc(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
       bufsize(:)=3*npw_block(:);bufdisp(:)=3*npw_disp(:)
       call xmpi_allgatherv(kg_k_loc,3*npw_k,kg_k,bufsize,bufdisp,spaceComm_band,ierr)
     else
       kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     end if
     do ilm=1,mpsang*mpsang
       ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       if (useylmgr>0) ylmgr_k(1:npw_k,1:3,ilm)=ylmgr(1+ikg:npw_k+ikg,1:3,ilm)
     end do

!    Compute (k+G) vectors
     ABI_ALLOCATE(kpg_k,(npw_nk,nkpg))
     if (nkpg>0) then
       call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_nk)
     end if
!    Allocate and compute the arrays phkxred and ph3d
     ABI_ALLOCATE(phkxred,(2,ncprj))
     do ia=iatom1,iatom2
       iatm=min(atindx_atm(ia),ncprj)
       arg=two_pi*(kpoint(1)*xred(1,ia)+kpoint(2)*xred(2,ia)+kpoint(3)*xred(3,ia))
       phkxred(1,iatm)=cos(arg);phkxred(2,iatm)=sin(arg)
     end do
     matblk=ncprj;if (nloalg(2)<=0) matblk=0
     ABI_ALLOCATE(ph3d,(2,npw_nk,matblk))
     if (matblk>0)then
!      Here, precomputation of ph3d
       if (npband_bandfft>1) then
!        Special treatment for band-fft //
         ABI_ALLOCATE(ph3d_tmp,(2,npw_k,matblk))
         call ph1d3d(1,ncprj,kg_k_loc,matblk,ncprj,npw_k,n1,n2,n3,phkxred,ph1d_atm,ph3d_tmp)
         ABI_ALLOCATE(ph3d_tmp_npw,(2,matblk,npw_k))
         ABI_ALLOCATE(ph3d_npw,(2,matblk,npw_nk))
         isize=2*matblk;bufsize(:)=isize*npw_block(:);bufdisp(:)=isize*npw_disp(:)
         do ipw=1,npw_k
           ph3d_tmp_npw(:,:,ipw)=ph3d_tmp(:,ipw,:)
         end do
         call xmpi_allgatherv(ph3d_tmp_npw,isize*npw_k,ph3d_npw,bufsize,bufdisp,spaceComm_band,ierr)
         do ipw=1,npw_nk
           ph3d(:,ipw,:)=ph3d_npw(:,:,ipw)
         end do
         ABI_DEALLOCATE(ph3d_npw)
         ABI_DEALLOCATE(ph3d_tmp_npw)
         ABI_DEALLOCATE(ph3d_tmp)
       else
         call ph1d3d(1,ncprj,kg_k,matblk,ncprj,npw_k,n1,n2,n3,phkxred,ph1d_atm,ph3d)
       end if
     else if (npband_bandfft>1) then
       MSG_ERROR('Band-fft parallelism +nloag(1)<0 forbidden !')
     end if

!    Compute nonlocal form factors ffnl at all (k+G)
     ABI_ALLOCATE(ffnl,(npw_nk,dimffnl,psps%lmnmax,ntypat0))
     if (npband_bandfft>1) then
!      Special treatment for band-fft //
       ABI_ALLOCATE(ffnl_tmp,(npw_k,dimffnl,psps%lmnmax,ntypat0))
       call mkffnl(psps%dimekb,dimffnl,ekb_atm,ffnl_tmp,ffspl_atm,&
&       gmet,gprimd,ider,idir0,indlmn_atm,kg_k_loc,kpg_k,kpoint,psps%lmnmax,&
&       psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,ntypat0,&
&       pspso_atm,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_k)
       ABI_ALLOCATE(ffnl_tmp_npw,(dimffnl,psps%lmnmax,ntypat0,npw_k))
       ABI_ALLOCATE(ffnl_npw,(dimffnl,psps%lmnmax,ntypat0,npw_nk))
       isize=dimffnl*psps%lmnmax*ntypat0
       bufsize(:)=isize*npw_block(:);bufdisp(:)=isize*npw_disp(:)
       do ipw=1,npw_k
         ffnl_tmp_npw(:,:,:,ipw)=ffnl_tmp(ipw,:,:,:)
       end do
       call xmpi_allgatherv(ffnl_tmp_npw,isize*npw_k,ffnl_npw,bufsize,bufdisp,spaceComm_band,ierr)
       do ipw=1,npw_nk
         ffnl(ipw,:,:,:)=ffnl_npw(:,:,:,ipw)
       end do
       ABI_DEALLOCATE(ffnl_npw)
       ABI_DEALLOCATE(ffnl_tmp_npw)
       ABI_DEALLOCATE(ffnl_tmp)
     else
       call mkffnl(psps%dimekb,dimffnl,ekb_atm,ffnl,ffspl_atm,&
&       gmet,gprimd,ider,idir0,indlmn_atm,kg_k,kpg_k,kpoint,psps%lmnmax,&
&       psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,ntypat0,&
&       pspso_atm,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_k)
     end if

!    No more need of kg_g_tmp
     if (npband_bandfft>1)  then
       ABI_DEALLOCATE(kg_k_loc)
     end if

!    Allocate arrays for a wave-function (or a block of WFs)
     ABI_ALLOCATE(cwavef,(2,npw_nk*my_nspinor*cg_bandpp))
     if (npband_bandfft>1) then
       isize=2*my_nspinor*cg_bandpp;bufsize(:)=isize*npw_block(:);bufdisp(:)=isize*npw_disp(:)
       isize=2*my_nspinor*npw_k*cg_bandpp;bufsize_wf(:)=isize
       do ii=1,npband_bandfft*cg_bandpp
         bufdisp_wf(ii)=(ii-1)*isize
       end do
     end if

!    Loop over bands or blocks of bands
     icgb=icg ; ibgb=ibg ; iband_start=1
     blocksz=npband_bandfft*cg_bandpp
     nblockbd=nband_k/blocksz
     nband_cprj_k=merge(nband_k/npband_bandfft,nband_k,cprj_band_distributed)
     do iblockbd=1,nblockbd
       iband_min=1+(iblockbd-1)*blocksz
       iband_max=iblockbd*blocksz

       if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband_min,iband_max,isppol,me_distrb)) then
         if (.not.cg_band_distributed) icgb=icgb+npw_k*my_nspinor*blocksz
         if (.not.cprj_band_distributed) ibgb=ibgb+my_nspinor*blocksz
         cycle
       end if

!      Extract wavefunction information
!      Special treatment for band-fft parallelism
       if (npband_bandfft>1) then
         !Transpose WF to get them in "FFT" representation
         ABI_ALLOCATE(cwavef_tmp,(2,npw_k*my_nspinor*blocksz))
         do ig=1,npw_k*my_nspinor*blocksz
           cwavef_tmp(1,ig)=cg(1,ig+icgb)
           cwavef_tmp(2,ig)=cg(2,ig+icgb)
         end do
         call xmpi_alltoallv(cwavef_tmp,bufsize_wf,bufdisp_wf,cwavef,bufsize,bufdisp,spaceComm_band,ierr)
         ABI_DEALLOCATE(cwavef_tmp)
         !Reorder WF according to cg_bandpp and/or spinor
         if (cg_bandpp>1.or.my_nspinor>1) then
           ABI_ALLOCATE(cwavef_tmp,(2,npw_nk*my_nspinor*blocksz))
           do ig=1,npw_nk*my_nspinor*blocksz
             cwavef_tmp(:,ig)=cwavef(:,ig)
           end do
           shift1=0
           do iwf2=1,cg_bandpp
             do ig=1,my_nspinor
               shift2=0
               do iwf1=1,npband_bandfft
                 npws=npw_block(iwf1)
                 ipw=shift2+(iwf2-1)*my_nspinor*npws+(ig-1)*npws
                 cwavef(:,shift1+1:shift1+npws)=cwavef_tmp(:,ipw+1:ipw+npws)
                 shift1=shift1+npws ; shift2=shift2+cg_bandpp*my_nspinor*npws
               end do
             end do
           end do
           ABI_DEALLOCATE(cwavef_tmp)
         end if
       else
         do ig=1,npw_k*my_nspinor*cg_bandpp
           cwavef(1,ig)=cg(1,ig+icgb)
           cwavef(2,ig)=cg(2,ig+icgb)
         end do
       end if

!      Compute scalar product of wavefunction with all NL projectors
       do ibp=1,cg_bandpp   ! Note: we suppose cp_bandpp=cprj_bandpp
         iwf1=1+(ibp-1)*npw_nk*my_nspinor;iwf2=ibp*npw_nk*my_nspinor
         icp1=1+(ibp-1)*my_nspinor;icp2=ibp*my_nspinor
         do jdir=istart,iend
           call getcprj(choice,cpopt,cwavef(:,iwf1:iwf2),cwaveprj(:,icp1:icp2),&
&           ffnl,jdir,indlmn_atm,istwf_k,kg_k,kpg_k,kpoint,psps%lmnmax,&
&           mgfft,mpi_enreg,ncprj,nattyp_atm,ngfft,nloalg,&
&           npw_nk,my_nspinor,ntypat0,phkxred,ph1d_atm,ph3d,ucvol,psps%useylm)
         end do
       end do
!      Export cwaveprj to big array cprj
       call pawcprj_put(atindx_atm,cwaveprj,cprj,ncprj,iband_start,ibgb,ikpt,iorder_cprj,isppol,&
&       mband_cprj,mkmem,natom,cprj_bandpp,nband_cprj_k,dimlmn,my_nspinor,nsppol,uncp,&
&       mpi_comm_band=spaceComm_band,to_be_gathered=(cg_band_distributed.and.(.not.cprj_band_distributed)))

       iband_start=iband_start+merge(cg_bandpp,blocksz,cprj_band_distributed)

!      End loop over bands
       icgb=icgb+npw_k*my_nspinor*blocksz
     end do

!    Shift array memory (if mkmem/=0)
     if (mkmem/=0) then
       ibg=ibg+my_nspinor*nband_cprj_k
       icg=icg+my_nspinor*nband_k*npw_k
       ikg=ikg+npw_k
     end if

!    End big k point loop
     ABI_DEALLOCATE(ffnl)
     ABI_DEALLOCATE(ph3d)
     ABI_DEALLOCATE(phkxred)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylmgr_k)
     ABI_DEALLOCATE(cwavef)
   end do
!  End loop over spins
 end do

 if ((iatom<=0).and.(choice==23)) then
   do iatom1=1,ncprj
     do ii=1,mcprj
       nlmn=cprj(iatom1,ii)%nlmn
       do kk=1,nlmn
         work(1:6)=cprj(iatom1,ii)%dcp(1,1:6,kk)
         call strconv(work,gprimd,work)
         cprj(iatom1,ii)%dcp(1,1:6,kk)=work(1:6)
         work(1:6)=cprj(iatom1,ii)%dcp(2,1:6,kk)
         call strconv(work,gprimd,work)
         cprj(iatom1,ii)%dcp(2,1:6,kk)=work(1:6)
       end do
     end do
   end do
 end if

!If needed, gather computed scalars
 if (.not.cg_band_distributed) then
   call pawcprj_mpi_sum(cprj,spaceComm_band,ierr)
 end if

!Deallocate temporary storage
 if (one_atom)  then
   ABI_DEALLOCATE(atindx_atm)
   ABI_DEALLOCATE(nattyp_atm)
   ABI_DEALLOCATE(ph1d_atm)
   ABI_DEALLOCATE(ekb_atm)
   ABI_DEALLOCATE(indlmn_atm)
   ABI_DEALLOCATE(ffspl_atm)
   ABI_DEALLOCATE(pspso_atm)
 end if
 nullify(atindx_atm,nattyp_atm,ph1d_atm,ekb_atm,indlmn_atm,ffspl_atm,pspso_atm)
 call pawcprj_free(cwaveprj)
 ABI_DATATYPE_DEALLOCATE(cwaveprj)
 ABI_DEALLOCATE(dimlmn)
 if (npband_bandfft>1) then
   ABI_DEALLOCATE(npw_block)
   ABI_DEALLOCATE(npw_disp)
   ABI_DEALLOCATE(bufsize)
   ABI_DEALLOCATE(bufdisp)
   ABI_DEALLOCATE(bufsize_wf)
   ABI_DEALLOCATE(bufdisp_wf)
 end if

 DBG_EXIT('COLL')

 end subroutine ctocprj
!!***

end module m_cgprj
!!***
