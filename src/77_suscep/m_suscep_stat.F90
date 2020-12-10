!!****m* ABINIT/m_suscep_stat
!! NAME
!! m_suscep_stat
!!
!! FUNCTION
!! Compute the susceptibility matrix
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (XG, AR, MB)
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

MODULE m_suscep_stat

 use defs_basis
 use m_xmpi
 use m_errors
 use m_abicore
 use m_distribfft

 use defs_abitypes, only : MPI_type
 use m_time,    only : timab
 use m_pawang,  only : pawang_type
 use m_pawtab,  only : pawtab_type
 use m_pawcprj, only : pawcprj_type, pawcprj_alloc, pawcprj_get, pawcprj_mpi_allgather, pawcprj_free
 use m_mpinfo,  only : destroy_mpi_enreg, initmpi_seq, proc_distrb_cycle
 use m_kg,      only : ph1d3d
 use m_gsphere, only : symg
 use m_fftcore, only : sphereboundary
 use m_fft,     only : fftpac, fourwf
 use m_spacepar,     only : symrhg
 use m_paw_finegrid, only : pawgylmg
 use m_paw_nhat,     only : pawsushat

 implicit none

 private
!!***

 public :: suscep_stat   ! Compute the susceptibility matrix

CONTAINS  !====================================================================================================
!!***

!!****f* m_suscep_stat/suscep_stat
!! NAME
!! suscep_stat
!!
!! FUNCTION
!! Compute the susceptibility matrix
!! from input wavefunctions, band occupations, and k point wts.
!! Include the usual sum-over-state terms, but also the
!! corrections due to the change of the Fermi level in the metallic
!! case, as well as implicit sum over higher lying conduction
!! states, thanks to the closure relation (referred to as an extrapolation).
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cg(2,mcg)=wf in G space
!!  cprj(natom,mcg*usecprj)= wave functions projected with non-local projectors:
!!                           cprj_nk(i)=<p_i|Cnk> where p_i is a non-local projector.
!!  dielar(7)=input parameters for dielectric matrix and susceptibility:
!!              diecut,dielng,diemac,diemix,diegap,dielam,diemixmag.
!!  dimcprj(natom*usepaw)=array of dimensions of array cprj (ordered by atom-type)
!!  doccde(mband*nkpt*nsppol)=derivative of occupancies wrt
!!           the energy for each band and k point
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  gbound_diel(2*mgfftdiel+8,2)=G sphere boundary for the dielectric matrix
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  irrzondiel(nfftdiel**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  lmax_diel=1+max. value of l angular momentum used for dielectric matrix
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mgfftdiel=maximum size of 1D FFTs, for the computation of the dielectric matrix
!!  mkmem=number of k points treated by this node
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum allowed value for npw
!!  natom=number of atoms in cell
!!  nband(nkpt*nsppol)=number of bands to be included in summation
!!   at each k point for each spin channel
!!  neglect_pawhat=1 if PAW contribution from hat density (compensation charge)
!!                 has to be neglected (to be used when only an estimation of
!!                 suscep. matrix has to be evaluated, i.e. for SCF precondictioning)
!!  nfftdiel=number of fft grid points for the computation of the diel matrix
!!  ngfftdiel(18)=contain all needed information about 3D FFT, for dielectric matrix,
!!    see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkpt=number of k points
!!  npwarr(nkpt)=number of planewaves and boundary planewaves
!!   at each k, for going from the WF sphere to the medium size FFT grid.
!!  npwdiel=third and fifth dimension of the susmat array.
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of symmetry elements in group (at least 1 for identity)
!!  ntypat=number of types of atoms in unit cell.
!!  occ(mband*nkpt*nsppol)=
!!          occupation numbers for each band (usually 2.0) at each k point
!!  occopt=option for occupancies
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  phnonsdiel(2,nfftdiel**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  ph1ddiel(2,3*(2*mgfftdiel+1)*natom*usepaw)=one-dimensional structure factor information
!!                                             for the dielectric matrix
!!  rprimd(3,3)=dimensional real space primitive translations
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrel(3,3,nsym)=symmetry matrices in real space (integers)
!!  tnons(3,nsym)=reduced nonsymmorphic translations
!!     (symrel and tnons are in terms of real space primitive translations)
!!  typat(natom)=type (integer) for each atom
!!  ucvol=unit cell volume (Bohr**3)
!!  unpaw=unit number for cprj PAW data (if used)
!!  usecprj= 1 if cprj array is stored in memory
!!  usepaw=flag for PAW
!!  usetimerev=1 if Time-Reversal symmetry has to be used when symmetrizing susceptibility
!!  wtk(nkpt)=k point weights (they sum to 1.0)
!!  ylmdiel(npwdiel,lmax_diel**2)= real spherical harmonics for each G and k point
!!                                 for the dielectric matrix
!!
!! OUTPUT
!!  susmat(2,npwdiel,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!
!! NOTES
!!  Case of non-collinear magnetism:
!!   In principle, we should compute 16 susceptibility matrix: chi0-(s1,s2),(s3,s4)
!!   (where s1, s2, s3,and s4 are spin indexes)...
!!   But, for the time being, the susceptibility is only used to compute the
!!   dielectric matrix within RPA approximation; in this approximation, only
!!   four susceptibilities are non-zero: chi0-(s1,s1),(s3,s3).
!!   They are stored in susmat(:,ipw1,1:2,ipw2,1:2)
!!
!! PARENTS
!!      m_vtorho
!!
!! CHILDREN
!!      destroy_mpi_enreg,fourwf,init_distribfft_seq,initmpi_seq,pawsushat
!!      timab
!!
!! SOURCE

subroutine suscep_stat(atindx,atindx1,cg,cprj,dielar,dimcprj,doccde,&
&  eigen,gbound_diel,gprimd,irrzondiel,istwfk,kg,&
&  kg_diel,lmax_diel,&
&  mband,mcg,mcprj,mgfftdiel,mkmem,mpi_enreg,mpw,natom,nband,&
&  neglect_pawhat,nfftdiel,ngfftdiel,nkpt,npwarr,&
&  npwdiel,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&
&  pawang,pawtab,phnonsdiel,ph1ddiel,rprimd,&
&  susmat,symafm,symrel,tnons,typat,ucvol,unpaw,usecprj,usepaw,usetimerev,&
&  wtk,ylmdiel)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmax_diel,mband,mcg,mcprj,mgfftdiel,mkmem,mpw,natom,neglect_pawhat
 integer,intent(in) :: nfftdiel,nkpt,npwdiel,nspden,nspinor,nsppol,nsym,ntypat,occopt
 integer,intent(in) :: unpaw,usecprj,usepaw,usetimerev
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: atindx(natom),atindx1(natom),dimcprj(natom*usepaw)
 integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
!no_abirules
!nfftdiel**(1-1/nsym) is 1 if nsym==1, and nfftdiel otherwise
 integer,intent(in) :: irrzondiel(nfftdiel**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4)),&
 & istwfk(nkpt)
 integer,intent(in) :: kg(3,mpw*mkmem),kg_diel(3,npwdiel),&
 & nband(nkpt*nsppol),ngfftdiel(18)
 integer,intent(in) :: npwarr(nkpt),symafm(nsym),symrel(3,3,nsym),typat(ntypat)
 real(dp),intent(in) :: cg(2,mcg),dielar(7)
 real(dp),intent(in) :: doccde(mband*nkpt*nsppol),eigen(mband*nkpt*nsppol)
 real(dp),intent(in) :: gprimd(3,3),occ(mband*nkpt*nsppol)
!nfftdiel**(1-1/nsym) is 1 if nsym==1, and nfftdiel otherwise
 real(dp),intent(in) :: phnonsdiel(2,nfftdiel**(1-1/nsym),(nspden/nsppol)-3*(nspden/4)),&
 &                                 tnons(3,nsym),wtk(nkpt)
 real(dp),intent(in) :: ph1ddiel(2,(3*(2*mgfftdiel+1)*natom)*usepaw),rprimd(3,3)
 real(dp),intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
 real(dp),intent(out) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 type(pawcprj_type) :: cprj(natom,mcprj*usecprj)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: bdtot_index,diag,extrap,i1,i2,i3,iband,ibg,icg,ier,ierr
 integer :: ifft,ii,ikg,ikpt,indx,iorder_cprj,ipw1,ipw2,isp,isp1,isp2
 integer :: ispinor,istwf_k,isym,j1,j2,j3,jj,jsp,k1,k2,k3
 integer :: my_nspinor,nband_k,nband_loc,ndiel1,ndiel2,ndiel3,ndiel4,ndiel5,ndiel6
 integer :: nkpg_diel,npw_k,npwsp,nspden_eff,nspden_tmp,nsym1,nsym2
 integer :: spaceComm,t1,t2,testocc
 real(dp) :: ai,ai2,ar,ar2,diegap,dielam,emax,invnsym
 real(dp) :: invnsym1,invnsym2,phi1,phi12,phi2,phr1,phr12
 real(dp) :: phr2,sumdocc,weight
 logical :: antiferro
 character(len=500) :: message
 type(MPI_type) :: mpi_enreg_diel

!arrays
 integer,allocatable :: gbound(:,:),kg_k(:,:),sym_g(:,:)
 integer,allocatable :: tmrev_g(:)
 real(dp) :: kpt_diel(3,1),tsec(2)
 real(dp),allocatable :: drhode(:,:,:),drhode_wk(:,:,:)
 real(dp),allocatable :: eig_diel(:),gylmg_diel(:,:,:),kpg_dum(:,:)
 real(dp),allocatable :: occ_deavg(:),ph3d_diel(:,:,:),phdiel(:,:,:)
 real(dp),allocatable :: phkxred_diel(:,:),rhoextrap(:,:,:,:),rhoextrg(:,:)
 real(dp),allocatable :: rhoextrr(:,:),sush(:),sussum(:),susvec(:,:,:)
 real(dp),allocatable :: suswk(:,:,:,:),zhpev1(:,:),zhpev2(:)
 type(pawcprj_type),allocatable :: cprj_k(:,:),cprj_loc(:,:)

! *************************************************************************

 call timab(740,1,tsec)
 call timab(741,1,tsec)

 ABI_CHECK(mkmem/=0,"mkmem==0 not supported anymore!")


!----- Initialisations -----------------------------------------------------------
!---------------------------------------------------------------------------------

 if (usecprj==0.and.usepaw==1) then
   write (message,'(3a)')&
&   ' cprj datastructure must be allocated !',ch10,&
&   ' Action: change pawusecp input keyword.'
   ABI_ERROR(message)
 end if

 if (mpi_enreg%paral_spinor==1) then
   message = ' not yet allowed for parallelization over spinors !'
   ABI_ERROR(message)
 end if

!Init mpicomm
 if(mpi_enreg%paral_kgb==1) then
   spaceComm=mpi_enreg%comm_kpt
 else
   spaceComm=mpi_enreg%comm_cell
 end if

!The dielectric stuff is performed in sequential mode.
!Set mpi_enreg_diel accordingly
 call initmpi_seq(mpi_enreg_diel)
 call init_distribfft_seq(MPI_enreg_diel%distribfft,'c',ngfftdiel(2),ngfftdiel(3),'all')

!testocc to be taken away
 testocc=1

 my_nspinor=max(1,nspinor/mpi_enreg%nproc_spinor)

 iorder_cprj=0 ! order for the cprj reading...

!Initialize some scalar quantities
 antiferro=(nsppol==1.and.nspden==2)
 nspden_eff=min(max(nsppol,nspden),2) ! Size for the computed part of susmat
 bdtot_index=0 ; icg=0 ; ibg=0
 ndiel1=ngfftdiel(1) ; ndiel2=ngfftdiel(2) ; ndiel3=ngfftdiel(3)

!ndiel4,ndiel5,ndiel6 are FFT dimensions, modified to avoid cache trashing
 ndiel4=ngfftdiel(4) ; ndiel5=ngfftdiel(5) ; ndiel6=ngfftdiel(6)
 diegap=dielar(5) ; dielam=dielar(6)
 extrap=0

!If dielam is too small, there is no extrapolation.
 if(dielam>1.0d-6)extrap=1

!Some stuff for symmetries
 nsym1=sum(symafm,mask=symafm==1)
 nsym2=nsym-nsym1
 invnsym =one/dble(nsym)
 invnsym1=one/dble(nsym1)
 invnsym2=one
!FIXME: make sure this is consistent with following code
!div by 0 for several v5 tests
 if (nsym2 > 0) invnsym2=one/dble(nsym2)

!Allocations
 ABI_MALLOC(occ_deavg,(mband))
 if(occopt>=3) then
   ABI_MALLOC(drhode,(2,npwdiel,nspden_eff))
 else
   ABI_MALLOC(drhode,(0,0,0))
 end if
 if(extrap==1) then
   ABI_MALLOC(rhoextrap,(ndiel4,ndiel5,ndiel6,nspinor))
 else
   ABI_MALLOC(rhoextrap,(0,0,0,0))
 end if

!zero the susceptibility matrix and other needed quantities
 susmat(:,:,:,:,:)=zero
 if(occopt>=3)then
   drhode(:,:,:)=zero
   sumdocc=zero
 end if

!PAW additional initializations
 if (usepaw==1) then
   ABI_MALLOC(gylmg_diel,(npwdiel,lmax_diel**2,ntypat))
   ABI_MALLOC(ph3d_diel,(2,npwdiel,natom))
   if (neglect_pawhat==0) then
     ABI_MALLOC(phkxred_diel,(2,natom))
     ABI_MALLOC(kpg_dum,(0,0))
     kpt_diel(1:3,1)=zero;phkxred_diel(1,:)=one;phkxred_diel(2,:)=zero;nkpg_diel=0
!    write(std_out,*) ' lmax_diel ', lmax_diel
     call pawgylmg(gprimd,gylmg_diel,kg_diel,kpg_dum,kpt_diel,lmax_diel,nkpg_diel,npwdiel,ntypat,pawtab,ylmdiel)
     call ph1d3d(1,natom,kg_diel,natom,natom,npwdiel,ndiel1,ndiel2,ndiel3,phkxred_diel,ph1ddiel,ph3d_diel)
     ABI_FREE(phkxred_diel)
     ABI_FREE(kpg_dum)
   else
     gylmg_diel=zero;ph3d_diel=one
   end if
 else
   ABI_MALLOC(gylmg_diel,(0,0,0))
   ABI_MALLOC(ph3d_diel,(0,0,0))
 end if

 call timab(741,2,tsec)



!--BIG loop over spins ------------------------------------------------------------
!---------------------------------------------------------------------------------

 do isp=1,nsppol
   ikg=0

   if(extrap==1)rhoextrap(:,:,:,:)=zero

!  --BIG loop over k-points --------------------------------------------------------
!  ---------------------------------------------------------------------------------

   do ikpt=1,nkpt

     nband_k=nband(ikpt+(isp-1)*nkpt)
     istwf_k=istwfk(ikpt)
     npw_k=npwarr(ikpt)

     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isp,mpi_enreg%me_kpt)) then
       bdtot_index=bdtot_index+nband_k
       cycle
     end if

     call timab(742,1,tsec)


     ABI_MALLOC(gbound,(2*mgfftdiel+8,2))
     ABI_MALLOC(kg_k,(3,npw_k))

     if (usepaw==1) then
       ABI_MALLOC(cprj_k,(natom,my_nspinor*nband_k))
       if (neglect_pawhat==0) then
         call pawcprj_alloc(cprj_k,0,dimcprj)
         if (mpi_enreg%nproc_band==1) then
           call pawcprj_get(atindx1,cprj_k,cprj,natom,1,ibg,ikpt,iorder_cprj,isp,&
&           mband,mkmem,natom,nband_k,nband_k,my_nspinor,nsppol,unpaw,&
&           mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         else
           nband_loc=nband_k/mpi_enreg%nproc_band
           ABI_MALLOC(cprj_loc,(natom,my_nspinor*nband_loc))
           call pawcprj_alloc(cprj_loc,0,dimcprj)
           call pawcprj_get(atindx1,cprj_loc,cprj,natom,1,ibg,ikpt,iorder_cprj,isp,&
&           mband/mpi_enreg%nproc_band,mkmem,natom,nband_loc,nband_loc,my_nspinor,nsppol,unpaw,&
&           mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
           call pawcprj_mpi_allgather(cprj_loc,cprj_k,natom,my_nspinor*nband_loc,mpi_enreg%bandpp,&
&           dimcprj,0,mpi_enreg%nproc_band,mpi_enreg%comm_band,ierr,rank_ordered=.true.)
           call pawcprj_free(cprj_loc)
           ABI_FREE(cprj_loc)
         end if
       else
         !call pawcprj_nullify(cprj_k)
       end if
     else
       ABI_MALLOC(cprj_k,(0,0))
     end if

     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     call sphereboundary(gbound,istwf_k,kg_k,mgfftdiel,npw_k)

     if(extrap==1)then
!      Compute inverse of average dielectric gap for each band
!      and multiply by occupation factor
       emax=maxval(eigen(1+bdtot_index:nband_k+bdtot_index))
       do iband=1,nband_k
         occ_deavg(iband)= occ(iband+bdtot_index)*dielam &
&         / ( emax-eigen(iband+bdtot_index)  + diegap )
       end do
     else
       occ_deavg(:)=zero
     end if

     call timab(742,2,tsec)
     call timab(743,1,tsec)

!    Compute the contribution of each k-point to susmat, rhoextrap, drhode and sumdocc.
     if(mpi_enreg%paral_kgb==1)then !Only this version is in parallel
!      Use either the simpler implementation
!      Should provide a test !!!
       call susk(atindx,bdtot_index,cg,cprj_k,doccde,drhode,eigen,extrap,gbound,&
&       gbound_diel,gylmg_diel,icg,ikpt,&
&       isp,istwf_k,kg_diel,kg_k,lmax_diel,mband,mcg,mgfftdiel,mpi_enreg,&
&       natom,nband_k,ndiel4,ndiel5,ndiel6,neglect_pawhat,ngfftdiel,nkpt,&
&       npwdiel,npw_k,nspden,nspden_eff,nspinor,nsppol,ntypat,occ,occopt,occ_deavg,&
&       pawang,pawtab,ph3d_diel,rhoextrap,sumdocc,&
&       susmat,typat,ucvol,usepaw,wtk)
     else
!      Or the more sophisticated one, needed to save memory.
       call suskmm(atindx,bdtot_index,cg,cprj_k,doccde,drhode,eigen,extrap,gbound,&
&       gbound_diel,gylmg_diel,icg,ikpt,&
&       isp,istwf_k,kg_diel,kg_k,lmax_diel,mband,mcg,mgfftdiel,mpi_enreg,&
&       natom,nband_k,ndiel4,ndiel5,ndiel6,neglect_pawhat,ngfftdiel,nkpt,&
&       npwdiel,npw_k,nspden,nspden_eff,nspinor,nsppol,ntypat,occ,occopt,occ_deavg,&
&       pawang,pawtab,ph3d_diel,rhoextrap,sumdocc,&
&       susmat,typat,ucvol,usepaw,wtk)
     end if

     call timab(743,2,tsec)

     ABI_FREE(gbound)
     ABI_FREE(kg_k)

     bdtot_index=bdtot_index+nband_k

     if (mkmem/=0) then
       ibg=ibg+my_nspinor*nband_k
       icg=icg+my_nspinor*npw_k*nband_k
       ikg=ikg+npw_k
     end if
     if (usepaw==1) then
       if (neglect_pawhat==0) then
         call pawcprj_free(cprj_k)
       end if
     end if
     ABI_FREE(cprj_k)

!    End loop on ikpt:  --------------------------------------------------------
   end do

!  Here include the contribution from the extrapolation to susmat,
!  diagonal part
   if(extrap==1)then

     call timab(744,1,tsec)

!    Transfer extrapolating density on augmented fft grid to
!    normal fft grid in real space.
!    Warning1 : if collinear magnetism, must treat only one spin at a time
!    Warning2 : if non-collinear magnetism, treat both spins
!    Warning3 : this is subtle for antiferro magnetism
     nspden_tmp=1;if (antiferro) nspden_tmp=2
     ABI_MALLOC(rhoextrr,(nfftdiel,nspden_tmp))
     ABI_MALLOC(rhoextrg,(2,nfftdiel))
     if (nspden==1.and.nspinor==2) rhoextrap(:,:,:,1)=rhoextrap(:,:,:,1)+rhoextrap(:,:,:,2)

     do ispinor=1,min(nspinor,nspden)
       jsp=isp+ispinor-1

       call fftpac(1,mpi_enreg_diel,1,ndiel1,ndiel2,ndiel3,ndiel4,ndiel5,ndiel6,&
&       ngfftdiel,rhoextrr(:,1),rhoextrap(:,:,:,ispinor),1)

!      Generate the density in reciprocal space, and symmetrize it
!      (note symrhg also make the reverse FFT, to get symmetrized density;
!      this is useless here, and should be made an option)
       call symrhg(1,gprimd,irrzondiel,mpi_enreg_diel,nfftdiel,nfftdiel,ngfftdiel,&
&       nspden_tmp,1,nsym,phnonsdiel,rhoextrg,rhoextrr,rprimd,symafm,symrel,tnons)

       do ipw2=1,npwdiel
         j1=kg_diel(1,ipw2) ; j2=kg_diel(2,ipw2) ; j3=kg_diel(3,ipw2)
!        static:    Only fills lower half of the matrix (here, the susceptibility matrix)
!        dynamical: fill all, will not affect susopt==2 for which extrap==0
         do ipw1=1,npwdiel
           i1=kg_diel(1,ipw1) ; i2=kg_diel(2,ipw1) ; i3=kg_diel(3,ipw1)
!          NOTE that there is a FFT folding (superposition) bias here
!          Should use kgindex, in the same spirit as in prcref
           k1=i1-j1; k1=modulo(k1,ndiel1)
           k2=i2-j2; k2=modulo(k2,ndiel2)
           k3=i3-j3; k3=modulo(k3,ndiel3)
           ifft=k1+1+ndiel1*(k2+ndiel2*k3)
           susmat(1,ipw1,jsp,ipw2,jsp)=susmat(1,ipw1,jsp,ipw2,jsp)+rhoextrg(1,ifft)
           susmat(2,ipw1,jsp,ipw2,jsp)=susmat(2,ipw1,jsp,ipw2,jsp)+rhoextrg(2,ifft)
         end do
       end do

     end do
     ABI_FREE(rhoextrg)
     ABI_FREE(rhoextrr)

     call timab(744,2,tsec)

   end if

!  End loop over spins ---------------------------------------------------------
 end do

 ABI_FREE(occ_deavg)
 ABI_FREE(rhoextrap)
 ABI_FREE(gylmg_diel)
 ABI_FREE(ph3d_diel)
 !end if

 call destroy_mpi_enreg(mpi_enreg_diel)

!-- Stuff for parallelism --------------------------------------------------------
!---------------------------------------------------------------------------------

 if(xmpi_paral==1)then
   call timab(746,1,tsec)
   ABI_MALLOC(sussum,(2*npwdiel*nspden*npwdiel*nspden))
!  Recreate full susmat on all proc.
!  This should be coded more efficiently,
!  since half of the matrix is still empty, and
!  it is spin-diagonal.
   sussum(:)=reshape(susmat(:,:,:,:,:),(/2*npwdiel*nspden*npwdiel*nspden/))
   call xmpi_sum(sussum,spaceComm,ierr)
   susmat(:,:,:,:,:)=reshape(sussum(:),(/2,npwdiel,nspden,npwdiel,nspden/))
   ABI_FREE(sussum)
!  Recreate full drhode on all proc.
   if(occopt>=3 .and. testocc==1)then
     call xmpi_sum(drhode,spaceComm,ierr)
!    Should use only one mpi-allreduce call instead of the three
     call xmpi_sum(sumdocc,spaceComm,ierr)
   end if
   call timab(746,2,tsec)
 end if

!-- Apply spatial hermitian/symmetries on spin-diagonal susceptibility matrix ----
!---------------------------------------------------------------------------------

 call timab(747,1,tsec)

!If antiferro magnetism, has to divide (spin-diagonal) susceptibility by 2 (due to dble occupations)
 if (antiferro) then
   do ipw2=1,npwdiel
     do ipw1=ipw2,npwdiel
       susmat(:,ipw1,1,ipw2,1)=half*susmat(:,ipw1,1,ipw2,1)
     end do
   end do
 end if

!Generate upper half of the spin-diagonal matrix (still the susceptibility matrix)
 do isp=1,nspden_eff
   do ipw2=2,npwdiel
     do ipw1=1,ipw2-1
       susmat(1,ipw1,isp,ipw2,isp)= susmat(1,ipw2,isp,ipw1,isp)
       susmat(2,ipw1,isp,ipw2,isp)=-susmat(2,ipw2,isp,ipw1,isp)
     end do
   end do
 end do

!Compute symmetric of G-vectors and eventual phases
!(either time-reversal or spatial symmetries)
 ABI_MALLOC(tmrev_g,(npwdiel))
 ABI_MALLOC(sym_g,(npwdiel,nsym))
 ABI_MALLOC(phdiel,(2,npwdiel,nsym))
 call symg(kg_diel,npwdiel,nsym,phdiel,sym_g,symrel,tmrev_g,tnons)

!Impose spatial symmetries to the spin-diagonal susceptibility matrix
 ABI_MALLOC(suswk,(2,npwdiel,npwdiel,nspden_eff))
 do isp=1,nspden_eff
   suswk(:,:,:,isp)=susmat(:,:,isp,:,isp) ! Temporary storage
 end do

 do isp=1,nspden_eff
   jsp=min(3-isp,nsppol)
   do ipw2=1,npwdiel
     do ipw1=1,npwdiel
       ar=suswk(1,ipw1,ipw2,isp)
       ai=suswk(2,ipw1,ipw2,isp)
       ar2=zero;ai2=zero
       if(nsym>1)then
         do isym=2,nsym
           t1=sym_g(ipw1,isym) ; t2=sym_g(ipw2,isym)
!          Not all symmetries are non-symmorphic. Should save time here ...
           phr1=phdiel(1,ipw1,isym) ; phi1=phdiel(2,ipw1,isym)
           phr2=phdiel(1,ipw2,isym) ; phi2=phdiel(2,ipw2,isym)
           phr12= phr1*phr2+phi1*phi2 ; phi12=phi1*phr2-phr1*phi2
           if (symafm(isym)==1) then
             ar=ar+suswk(1,t1,t2,isp)*phr12-suswk(2,t1,t2,isp)*phi12
             ai=ai+suswk(2,t1,t2,isp)*phr12+suswk(1,t1,t2,isp)*phi12
           else
             ar2=ar2+suswk(1,t1,t2,jsp)*phr12-suswk(2,t1,t2,jsp)*phi12
             ai2=ai2+suswk(2,t1,t2,jsp)*phr12+suswk(1,t1,t2,jsp)*phi12
           end if
         end do
       end if
       if (antiferro) then
         susmat(1,ipw1,1,ipw2,1)=ar*invnsym1
         susmat(2,ipw1,1,ipw2,1)=ai*invnsym1
         susmat(1,ipw1,2,ipw2,2)=ar2*invnsym2
         susmat(2,ipw1,2,ipw2,2)=ai2*invnsym2
       else
         susmat(1,ipw1,isp,ipw2,isp)=(ar+ar2)*invnsym
         susmat(2,ipw1,isp,ipw2,isp)=(ai+ai2)*invnsym
       end if
     end do
   end do
 end do
 ABI_FREE(suswk)


!--  Add contribibution to susceptibility due to change of Fermi level -----------
!---------------------------------------------------------------------------------

 if (occopt>=3.and.testocc==1) then

!  Impose spatial symmetries to drhode
   ABI_MALLOC(drhode_wk,(2,npwdiel,nspden_eff))
   do isp=1,nspden_eff
     jsp=min(3-isp,nsppol)
     do ipw1=1,npwdiel
       ar=drhode(1,ipw1,isp)
       ai=drhode(2,ipw1,isp)
       ar2=zero;ai2=zero
       if (nsym>1) then
         do isym=2,nsym
           t1=sym_g(ipw1,isym)
!          Not all symmetries are non-symmorphic. Should save time here ...
           phr1=phdiel(1,ipw1,isym);phi1=phdiel(2,ipw1,isym)
           if (symafm(isym)==1) then
             ar=ar+drhode(1,t1,isp)*phr1-drhode(2,t1,isp)*phi1
             ai=ai+drhode(2,t1,isp)*phr1+drhode(1,t1,isp)*phi1
           else
             ar2=ar2+drhode(1,t1,jsp)*phr1-drhode(2,t1,jsp)*phi1
             ai2=ai2+drhode(2,t1,jsp)*phr1+drhode(1,t1,jsp)*phi1
           end if
         end do
       end if
       if (antiferro) then  ! 1/2 factor due to (dble) occupations
         drhode_wk(1,ipw1,1)=half*ar*invnsym1
         drhode_wk(2,ipw1,1)=half*ai*invnsym1
         drhode_wk(1,ipw1,2)=half*ar2*invnsym2
         drhode_wk(2,ipw1,2)=half*ai2*invnsym2
       else
         drhode_wk(1,ipw1,isp)=(ar+ar2)*invnsym
         drhode_wk(2,ipw1,isp)=(ai+ai2)*invnsym
       end if
     end do
   end do

!  Add contribution to non-diagonal susceptibility
!  Presently fills complete susceptibility matrix, not only lower half
   weight=one/sumdocc
   do isp2=1,nspden_eff
     do ipw2=1,npwdiel
       do isp1=1,nspden_eff
         do ipw1=1,npwdiel
           susmat(1,ipw1,isp1,ipw2,isp2)=susmat(1,ipw1,isp1,ipw2,isp2)- &
&           weight*( drhode_wk(1,ipw1,isp1)*drhode_wk(1,ipw2,isp2)  &
&           +drhode_wk(2,ipw1,isp1)*drhode_wk(2,ipw2,isp2) )
           susmat(2,ipw1,isp1,ipw2,isp2)=susmat(2,ipw1,isp1,ipw2,isp2)- &
&           weight*( drhode_wk(2,ipw1,isp1)*drhode_wk(1,ipw2,isp2)  &
&           -drhode_wk(1,ipw1,isp1)*drhode_wk(2,ipw2,isp2) )
         end do
       end do
     end do
   end do
   ABI_FREE(drhode_wk)

 end if
 !if (occopt>=3)  then
 ABI_FREE(drhode)
 !end if


!--- Impose the time-reversal symmetry to the susceptibility matrix --------------
!---------------------------------------------------------------------------------

 if (usetimerev==1) then
   ABI_MALLOC(suswk,(2,npwdiel,npwdiel,1))

!  Impose the time-reversal symmetry to the spin-diagonal susceptibility matrix
   do isp=1,nspden_eff
     suswk(:,:,:,1)=susmat(:,:,isp,:,isp) ! Temporary storage
     do ipw2=1,npwdiel
       t2=tmrev_g(ipw2)
       do ipw1=1,npwdiel
         t1=tmrev_g(ipw1)
         susmat(1,ipw1,isp,ipw2,isp)=half*(suswk(1,ipw1,ipw2,1)+suswk(1,t1,t2,1))
         susmat(2,ipw1,isp,ipw2,isp)=half*(suswk(2,ipw1,ipw2,1)-suswk(2,t1,t2,1))
       end do
     end do
   end do

!  Impose the time-reversal symmetry to the off-diagonal susceptibility matrix
   if (nspden_eff/=1.and.occopt>=3.and.testocc==1) then
     suswk(:,:,:,1)=susmat(:,:,1,:,2) ! Temporary storage
     do ipw2=1,npwdiel
       t2=tmrev_g(ipw2)
       do ipw1=1,npwdiel
         t1=tmrev_g(ipw1)
         ar=half*(suswk(1,ipw1,ipw2,1)+suswk(1,t1,t2,1))
         ai=half*(suswk(2,ipw1,ipw2,1)-suswk(2,t1,t2,1))
         susmat(1,ipw1,1,ipw2,2)= ar
         susmat(2,ipw1,1,ipw2,2)= ai
         susmat(1,ipw1,2,ipw2,1)= ar
         susmat(2,ipw1,2,ipw2,1)=-ai
       end do
     end do
   end if
   ABI_FREE(suswk)
 end if

 ABI_FREE(phdiel)
 ABI_FREE(sym_g)
 ABI_FREE(tmrev_g)


!-- The full susceptibility matrix is computed -----------------------------------
!-- Now, eventually diagonalize it and stop --------------------------------------
!---------------------------------------------------------------------------------

!Must turn on this flag to make the diagonalisation
 diag=0
 if(diag==1)then

   npwsp=npwdiel*nspden_eff
   ABI_MALLOC(sush,(npwsp*(npwsp+1)))
   ABI_MALLOC(susvec,(2,npwsp,npwsp))
   ABI_MALLOC(eig_diel,(npwsp))
   ABI_MALLOC(zhpev1,(2,2*npwsp-1))
   ABI_MALLOC(zhpev2,(3*npwsp-2))
   ier=0

!  Store the susceptibility matrix in proper mode before calling zhpev
   indx=1
   do ii=1,npwdiel
     do jj=1,ii
       sush(indx  )=susmat(1,jj,1,ii,1)
       sush(indx+1)=susmat(2,jj,1,ii,1)
       indx=indx+2
     end do
   end do

!  If spin-polarized, need to store other parts of the matrix
   if(nspden_eff/=1)then
     do ii=1,npwdiel
!      Here, spin-flip contribution
       do jj=1,npwdiel
         sush(indx  )=susmat(1,jj,1,ii,2)
         sush(indx+1)=susmat(2,jj,1,ii,2)
         indx=indx+2
       end do
!      Here spin down-spin down upper matrix
       do jj=1,ii
         sush(indx  )=susmat(1,jj,2,ii,2)
         sush(indx+1)=susmat(2,jj,2,ii,2)
         indx=indx+2
       end do
     end do
   end if

   call ZHPEV ('V','U',npwsp,sush,eig_diel,susvec,npwsp,zhpev1,&
&   zhpev2,ier)

   write(std_out,*)' suscep_stat : print eigenvalues of the susceptibility matrix'
   do ii=1,npwsp
     write(std_out,'(i5,es16.6)' )ii,eig_diel(ii)
   end do

   ABI_FREE(sush)
   ABI_FREE(susvec)
   ABI_FREE(eig_diel)
   ABI_FREE(zhpev1)
   ABI_FREE(zhpev2)
   ABI_ERROR("Stopping here!")
 end if

 call timab(747,2,tsec)
 call timab(740,2,tsec)

end subroutine suscep_stat
!!***

!!****f* m_suscep_stat/susk
!! NAME
!! susk
!!
!! FUNCTION
!! Compute the contribution of one k point to the susceptibility matrix
!! from input wavefunctions, band occupations, and k point wts.
!! Include the usual sum-over-state terms, but also the
!! corrections due to the change of the Fermi level in the metallic
!! case, as well as implicit sum over higher lying conduction
!! states, thanks to the closure relation (referred to as an extrapolation).
!! Compared to the routine suskmm, there is no particular attention
!! to the use of the memory, so the code is simpler.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms
!!  bdtot_index=index for the number of the band
!!  cg(2,mcg)=wfs in G space
!!  cprj_k(natom,nspinor*nband_k)= wave functions projected with non-local projectors:
!!                                 cprj_k=<p_i|Cnk> where p_i is a non-local projector.
!!  doccde(mband*nkpt*nsppol)=derivative of occupancies wrt
!!           the energy for each band and k point
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  extrap: if==1, the closure relation (an extrapolation) must be used
!!  gbound(2*mgfftdiel+8,2)=G sphere boundary for going from WF sphere to
!!      medium size FFT grid
!!  gbound_diel(2*mgfftdiel+8,2)=G sphere boundary for going from medium size
!!      FFT grid to small sphere.
!!  gylmg_diel(npwdiel,lmax_diel,ntypat*usepaw)= -PAW only- Fourier transform of g_l(r).Y_ml(r) shape functions
!!                                               for dielectric matrix
!!  icg=index for cg
!!  ikpt=number of the k point
!!  isp=number of the current spin
!!  istwf_k=input option parameter that describes the storage of wfs
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  kg_k(3,npw_k)=coordinates of planewaves in basis sphere.
!!  lmax_diel=1+max. value of l angular momentum used for dielectric matrix
!!  mband=maximum number of bands
!!  mcg=dimension of cg
!!  mgfftdiel=maximum size of 1D FFTs, for the computation of
!!     the dielectric matrix
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell
!!  nband_k=number of bands at this k point for that spin polarization
!!  ndiel4,ndiel5,ndiel6= FFT dimensions, modified to avoid cache trashing
!!  neglect_pawhat=1 if PAW contribution from hat density (compensation charge)
!!                 has to be neglected (to be used when only an estimation of
!!                 suscep. matrix has to be evaluated, i.e. for SCF precondictioning)
!!  ngfftdiel(18)=contain all needed information about 3D FFT, for dielectric matrix,
!!    see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkpt=number of k points
!!  npwdiel=third and fifth dimension of the susmat array.
!!  npw_k=number of plane waves at this k point
!!  nspden=number of spin-density components
!!  nspden_eff=number of spin-density components actually computed in sussceptibility
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  occ(mband*nkpt*nsppol)=
!!          occupation numbers for each band (usually 2.0) at each k point
!!  occopt=option for occupancies
!!  occ_deavg(mband)=factor for extrapolation (occup. divided by an energy gap)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ph3d_diel(2,npwdiel,natom*usepaw)=3-dim structure factors, for each atom and plane wave, for dielectric matrix
!!  typat(natom)=type (integer) for each atom
!!  ucvol=unit cell volume (Bohr**3)
!!  usepaw=flag for PAW
!!  wtk(nkpt)=k point weights (they sum to 1.0)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! These quantities are accumulated in this routine:
!!  drhode(2,npwdiel,nspden_eff)=weighted density, needed to compute the
!!   effect of change of fermi energy
!!  rhoextrap(ndiel4,ndiel5,ndiel6,nspinor)=density-like array, needed for the
!!   extrapolation procedure.
!!  sumdocc=sum of weighted occupation numbers, needed to compute the
!!   effect of change of fermi energy
!!  susmat(2,npwdiel,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!
!! NOTES
!! Band-fft parallel treatment: Each processor will treat his own band, but susmat will be known by all.
!! This means that cg will not have the same meaning in sequential or parallel mode.
!! In parallel mode, it will contain the set of all bands treated by the currrent processor.
!! To achieve this, the argument cg has been replaced by cg_mpi, with the "target" attribute.
!! In sequential mode, the pointer cg will point towards cg_mpi. In parallel mode, cg will point
!! to a new array cg_local, containing the bands treated by the currrent processor.
!! This allows to minimize the overhead incurred by the parallelization  of the sequential version.
!! A similar treatment is performed on kg_k, npw_k.
!! A future version might have objects like kg_k_gather as arguments, instead of computing them.
!! This is in slight violation of programming rules, but I think it is safe, since the pointers remain local
!! GZ
!!
!! PARENTS
!!      m_suscep_stat
!!
!! CHILDREN
!!      destroy_mpi_enreg,fourwf,init_distribfft_seq,initmpi_seq,pawsushat
!!      timab
!!
!! SOURCE

subroutine susk(atindx,bdtot_index,cg_mpi,cprj_k,doccde,drhode,eigen,extrap,gbound,&
&  gbound_diel,gylmg_diel,icg_mpi,ikpt,isp,istwf_k,kg_diel,kg_k_mpi,&
&  lmax_diel,mband,mcg,mgfftdiel,mpi_enreg,&
&  natom,nband_k,ndiel4,ndiel5,ndiel6,neglect_pawhat,ngfftdiel,nkpt,&
&  npwdiel,npw_k_mpi,nspden,nspden_eff,nspinor,nsppol,ntypat,occ,occopt,occ_deavg,&
&  pawang,pawtab,ph3d_diel,rhoextrap,sumdocc,&
&  susmat,typat,ucvol,usepaw,wtk)

!Arguments ------------------------------------
!This type is defined in defs_mpi
!scalars
 integer,intent(in) :: bdtot_index,extrap,ikpt,isp,istwf_k,lmax_diel,mband,mcg
 integer,intent(in) :: mgfftdiel,natom,nband_k,ndiel4,ndiel5,ndiel6,neglect_pawhat
 integer,intent(in) :: nkpt,npwdiel,nspden,nspden_eff,nspinor,nsppol
 integer,intent(in) :: ntypat,occopt,usepaw
 integer,intent(in),target :: icg_mpi,npw_k_mpi
 real(dp),intent(in) :: ucvol
 real(dp),intent(inout) :: sumdocc
 type(MPI_type),intent(in) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: atindx(natom),gbound_diel(2*mgfftdiel+8,2)
 integer,intent(in) :: kg_diel(3,npwdiel),ngfftdiel(18),typat(natom)
 integer,intent(in),target :: kg_k_mpi(3,npw_k_mpi)
 integer,intent(inout) :: gbound(2*mgfftdiel+8,2)
 integer,pointer :: kg_k(:,:)
 real(dp),intent(in) :: doccde(mband*nkpt*nsppol),eigen(mband*nkpt*nsppol)
 real(dp),intent(in) :: gylmg_diel(npwdiel,lmax_diel**2,ntypat*usepaw)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),occ_deavg(mband)
 real(dp),intent(in) :: ph3d_diel(2,npwdiel,natom*usepaw),wtk(nkpt)
 real(dp),intent(in),target :: cg_mpi(2,mcg)
 real(dp),intent(inout) :: drhode(2,npwdiel,nspden_eff)
 real(dp),intent(inout) :: rhoextrap(ndiel4,ndiel5,ndiel6,nspinor)
 real(dp),intent(inout) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 type(pawcprj_type) :: cprj_k(natom,nspinor*nband_k*usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
! real(dp), allocatable :: cg_disk(:,:)
!Local variables for MPI
!scalars
 integer :: blocksize,i1,i2,i3,iband,iband_loc,ibd1,ibd2,ibdblock,ier
 integer :: iproc,iproc_fft,ipw,ipw1,ipw2,isp1,isp2,ispinor,iwf,jsp,me_bandfft
 integer :: nbdblock,ndatarecv,ndiel1,ndiel2,ndiel3
 integer :: sizemax_per_proc,spaceComm,testocc,tim_fourwf
 integer,pointer :: icg,npw_k
 integer,target :: icg_loc=0,npw_k_loc,npw_tot
 real(dp) :: eigdiff,occdiff,tolocc,weight,wght1,wght2
 type(MPI_type) :: mpi_enreg_diel
!arrays
 integer,allocatable :: band_loc(:),kg_k_gather(:,:),npw_per_proc(:),rdispls(:)
 integer,allocatable :: rdispls_all(:),rdisplsloc(:),recvcounts(:)
 integer,allocatable :: recvcountsloc(:),sdispls(:),sdisplsloc(:),sendcounts(:)
 integer,allocatable :: sendcountsloc(:)
 integer,allocatable,target :: kg_k_gather_all(:,:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: cwavef(:,:),cwavef_alltoall(:,:)
 real(dp),allocatable :: cwavef_alltoall_gather(:,:),dummy(:,:),rhoaug(:,:,:)
 real(dp),allocatable :: susmat_mpi(:,:,:)
 real(dp),allocatable :: wfprod(:,:),wfraug(:,:,:,:),wfrspa(:,:,:,:,:,:)
 real(dp),allocatable,target :: cg_local(:,:)
 real(dp),pointer :: cg(:,:)
 logical,allocatable :: treat_band(:)

! *************************************************************************

!DEBUG
!write(std_out,*)' susk : enter '
!if(.true.)stop
!ENDDEBUG

 call timab(750,1,tsec)
 call timab(751,1,tsec)

 ndiel1=ngfftdiel(1) ; ndiel2=ngfftdiel(2) ; ndiel3=ngfftdiel(3)

!The dielectric stuff is performed in sequential mode.
!Set mpi_enreg_diel accordingly
 call initmpi_seq(mpi_enreg_diel)
 call init_distribfft_seq(mpi_enreg_diel%distribfft,'c',ngfftdiel(2),ngfftdiel(3),'all')
 me_bandfft=xmpi_comm_rank(mpi_enreg%comm_bandfft)

 testocc=1
!DEBUG
!write(std_out,*)' susk : set testocc to 0 '
!testocc=0
!write(std_out,*)' susk : set extrap to 0 '
!extrap=0
!ENDDEBUG

!Allocations, initializations
 ABI_MALLOC(rhoaug,(ndiel4,ndiel5,ndiel6))
 ABI_MALLOC(wfraug,(2,ndiel4,ndiel5,ndiel6))
 ABI_MALLOC(wfprod,(2,npwdiel))
 ABI_MALLOC(wfrspa,(2,ndiel4,ndiel5,ndiel6,nspinor,mband))
 ABI_MALLOC(dummy,(2,1))
 wfrspa(:,:,:,:,:,:)=zero
 ABI_MALLOC(treat_band,(nband_k))
 treat_band(:)=.true.
 isp1=isp;isp2=isp
 if (nspden_eff==2.and.nspinor==2) isp2=isp+1

!BAND-FFT parallelism
 if (mpi_enreg%paral_kgb==1) then
   treat_band(:)=.false.
!  We gather the wavefunctions treated by this proc in cg_local
   spaceComm=mpi_enreg%comm_band
   blocksize=mpi_enreg%nproc_band
   nbdblock=nband_k/blocksize
   ABI_MALLOC(sdispls,(blocksize))
   ABI_MALLOC(sdisplsloc,(blocksize))
   ABI_MALLOC(sendcounts,(blocksize))
   ABI_MALLOC(sendcountsloc,(blocksize))
   ABI_MALLOC(rdispls,(blocksize))
   ABI_MALLOC(rdisplsloc,(blocksize))
   ABI_MALLOC(recvcounts,(blocksize))
   ABI_MALLOC(recvcountsloc,(blocksize))
!  First gather the kg_k in kg_k_gather_all
   npw_k_loc=npw_k_mpi
   call xmpi_allgather(npw_k_loc,recvcounts,spaceComm,ier)
   rdispls(1)=0
   do iproc=2,blocksize
     rdispls(iproc)=rdispls(iproc-1)+recvcounts(iproc-1)
   end do
   ndatarecv=rdispls(blocksize)+recvcounts(blocksize)
   ABI_MALLOC(kg_k_gather,(3,ndatarecv))
   recvcountsloc(:)=recvcounts(:)*3
   rdisplsloc(:)=rdispls(:)*3
   call xmpi_allgatherv(kg_k_mpi,3*npw_k_loc,kg_k_gather,recvcountsloc(:),rdisplsloc,spaceComm,ier)
   ABI_MALLOC(npw_per_proc,(mpi_enreg%nproc_fft))
   ABI_MALLOC(rdispls_all,(mpi_enreg%nproc_fft))
   spaceComm=mpi_enreg%comm_fft
   call xmpi_allgather(ndatarecv,npw_per_proc,spaceComm,ier)
   rdispls_all(1)=0
   do iproc=2,mpi_enreg%nproc_fft
     rdispls_all(iproc)=rdispls_all(iproc-1)+npw_per_proc(iproc-1)
   end do
   npw_tot=rdispls_all(mpi_enreg%nproc_fft)+npw_per_proc(mpi_enreg%nproc_fft)
   ABI_MALLOC(kg_k_gather_all,(3,npw_tot))
   call xmpi_allgatherv(kg_k_gather,3*ndatarecv,kg_k_gather_all,3*npw_per_proc(:),3*rdispls_all,spaceComm,ier)
!  At this point kg_k_gather_all contains all the kg
   if(allocated(cwavef))  then
     ABI_FREE(cwavef)
   end if
   ABI_MALLOC(cwavef,(2,npw_k_loc*nspinor*blocksize))
   sizemax_per_proc=nband_k/(mpi_enreg%nproc_band*mpi_enreg%nproc_fft)+1
   ABI_MALLOC(band_loc,(sizemax_per_proc))
   ABI_MALLOC(cg_local,(2,sizemax_per_proc*npw_tot*nspinor))
   iband_loc=0
   do ibdblock=1,nbdblock
     cwavef(:,1:npw_k_loc*nspinor*blocksize)=&
&     cg_mpi(:,1+(ibdblock-1)*npw_k_loc*nspinor*blocksize+icg_mpi:ibdblock*npw_k_loc*nspinor*blocksize+icg_mpi)
     sendcounts(:)=npw_k_loc
     do iproc=1,blocksize
       sdispls(iproc)=(iproc-1)*npw_k_loc
     end do
     ABI_MALLOC(cwavef_alltoall,(2,ndatarecv*nspinor))
     recvcountsloc(:)=recvcounts(:)*2*nspinor
     rdisplsloc(:)=rdispls(:)*2*nspinor
     sendcountsloc(:)=sendcounts(:)*2*nspinor
     sdisplsloc(:)=sdispls(:)*2*nspinor
     call timab(547,1,tsec)
     spaceComm=mpi_enreg%comm_band
     call xmpi_alltoallv(cwavef,sendcountsloc,sdisplsloc,cwavef_alltoall,recvcountsloc,rdisplsloc,spaceComm,ier)
     call timab(547,2,tsec)
     ABI_MALLOC(cwavef_alltoall_gather,(2,npw_tot*nspinor))
     blocksize=mpi_enreg%nproc_band
     spaceComm=mpi_enreg%comm_fft
     call xmpi_allgatherv(cwavef_alltoall,2*nspinor*ndatarecv,cwavef_alltoall_gather,&
&     2*nspinor*npw_per_proc,2*nspinor*rdispls_all,spaceComm,ier)
     iproc_fft=modulo(ibdblock-1,mpi_enreg%nproc_fft)
     if(mpi_enreg%me_fft==iproc_fft) then !All nproc_band procs of index me_fft will treat these bands
       iband_loc=iband_loc+1
       iband=1+mpi_enreg%me_band+mpi_enreg%nproc_band*mpi_enreg%me_fft+(iband_loc-1)*mpi_enreg%nproc_fft*mpi_enreg%nproc_band
       treat_band(iband)=.true.
       band_loc(iband_loc)=iband
       cg_local(:,1+(iband_loc-1)*npw_tot*nspinor:iband_loc*npw_tot*nspinor)=cwavef_alltoall_gather(:,1:npw_tot*nspinor)
     end if
     ABI_FREE(cwavef_alltoall_gather)
     ABI_FREE(cwavef_alltoall)
   end do
!  On exit:
!  npw_tot will be npw
!  kg_k_gather_all will be kg_k
!  cg_local will be cg
!  icg will be zero
   npw_k=>npw_tot
   kg_k=>kg_k_gather_all(:,:)
   cg=>cg_local(:,:)
   icg=>icg_loc
   call sphereboundary(gbound,istwf_k,kg_k,mgfftdiel,npw_k)
   ABI_FREE(npw_per_proc)
   ABI_FREE(rdispls_all)
   ABI_FREE(sendcounts)
   ABI_FREE(recvcounts)
   ABI_FREE(sdispls)
   ABI_FREE(rdispls)
   ABI_FREE(sendcountsloc)
   ABI_FREE(sdisplsloc)
   ABI_FREE(recvcountsloc)
   ABI_FREE(rdisplsloc)
   ABI_FREE(kg_k_gather)
   ABI_FREE(cwavef)
!  Because they will be summed over all procs, and arrive on input, rescale drhode and rhoextrap
   if(occopt>=3)drhode(:,:,isp1:isp2)=drhode(:,:,isp1:isp2)/real(mpi_enreg%nproc_fft*mpi_enreg%nproc_band,dp)
   if(extrap==1)rhoextrap(:,:,:,:)=rhoextrap(:,:,:,:)/real(mpi_enreg%nproc_fft*mpi_enreg%nproc_band,dp)
   do i1=isp1,isp2
     susmat(:,:,i1,:,i1)=susmat(:,:,i1,:,i1)/real(mpi_enreg%nproc_fft*mpi_enreg%nproc_band,dp)
   end do

!  No BAND-FFT parallelism
 else ! use argument variables
   cg=>cg_mpi
   kg_k=>kg_k_mpi
   npw_k=>npw_k_mpi
   icg=>icg_mpi
 end if
 iband_loc=0

 call timab(751,2,tsec)
 call timab(752,1,tsec)

!Loop over bands to fft and store Fourier transform of wavefunction
 ABI_MALLOC(cwavef,(2,npw_k))
 do iband=1,nband_k
   if(.not. treat_band(iband))  cycle ! I am not treating this band (only for the parallel case)
   iband_loc=iband_loc+1

!  Loop on spinorial components
   do ispinor=1,nspinor
     iwf=(ispinor-1)*npw_k+(iband_loc-1)*npw_k*nspinor+icg
     jsp=isp+ispinor-1;if (nspden_eff==1) jsp=isp

!    Obtain Fourier transform in fft box
     tim_fourwf=8
     cwavef(:,1:npw_k)=cg(:,1+iwf:npw_k+iwf)
     call fourwf(1,rhoaug,cwavef,dummy,wfraug,gbound,gbound,&
&     istwf_k,kg_k,kg_k,mgfftdiel,mpi_enreg_diel,1,ngfftdiel,npw_k,1,ndiel4,ndiel5,ndiel6,&
&     0,tim_fourwf,weight,weight)

     wfrspa(:,:,:,:,ispinor,iband)=wfraug(:,:,:,:)

     if( (occopt>=3 .and. testocc==1) .or. extrap==1 )then
!      In the case of metallic occupation, or if the extrapolation
!      over higher bands is included, must compute the
!      Fourier transform of the density of each band, then
!      generate the part of the susceptibility matrix due
!      varying occupation numbers.

       weight=-two*occ_deavg(iband)*wtk(ikpt)/ucvol
       do i3=1,ndiel3
         do i2=1,ndiel2
           do i1=1,ndiel1
             wfraug(1,i1,i2,i3)=wfraug(1,i1,i2,i3)**2+wfraug(2,i1,i2,i3)**2
             wfraug(2,i1,i2,i3)=zero
           end do
         end do
!        If extrapolation, accumulate density in real space
         if(extrap==1.and.usepaw==0)then
           do i2=1,ndiel2
             do i1=1,ndiel1
               rhoextrap(i1,i2,i3,ispinor)=rhoextrap(i1,i2,i3,ispinor)+weight*wfraug(1,i1,i2,i3)
             end do
           end do
         end if
       end do

!      In case of PAW, add compensation charge contribution
       if (usepaw==1.and.extrap==1.and.neglect_pawhat==0) then
         call pawsushat(atindx,cprj_k,gbound_diel,gylmg_diel,iband,iband,ispinor,ispinor,1,kg_diel,&
&         lmax_diel,mgfftdiel,natom,nband_k,ndiel4,ndiel5,ndiel6,&
&         ngfftdiel,npwdiel,nspinor,ntypat,1,&
&         pawang,pawtab,ph3d_diel,typat,dummy,wfraug,&
&         mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
         rhoextrap(:,:,:,ispinor)=rhoextrap(:,:,:,ispinor)+weight*wfraug(1,:,:,:)
       end if

!      Performs the Fourier Transform of the density of the band,
!      and store it in wfprod
       tim_fourwf=9
       call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&       1,kg_diel,kg_diel,mgfftdiel,mpi_enreg_diel,1,ngfftdiel,1,npwdiel,&
&       ndiel4,ndiel5,ndiel6,3,tim_fourwf,weight,weight)
!      In case of PAW, add compensation charge contribution if not already done
       if (usepaw==1.and.extrap==0.and.neglect_pawhat==0) then
         call pawsushat(atindx,cprj_k,gbound_diel,gylmg_diel,ibd1,ibd2,ispinor,ispinor,1,kg_diel,&
&         lmax_diel,mgfftdiel,natom,nband_k,ndiel4,ndiel5,ndiel6,&
&         ngfftdiel,npwdiel,nspinor,ntypat,0,&
&         pawang,pawtab,ph3d_diel,typat,wfprod,dummy,&
&         mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
       end if

!      Perform now the summation of terms related to direct change of eigenvalues
!      or extrapolation over higher bands
       wght1=zero ; wght2=zero
       if(occopt>=3 .and. testocc==1) wght1=doccde(iband+bdtot_index)*wtk(ikpt)/ucvol
       if(extrap==1) wght2=two*occ_deavg(iband)*wtk(ikpt)/ucvol
       weight=wght1+wght2

       if (abs(weight)>tol12) then
         do ipw2=1,npwdiel
!          Only fills lower half of the matrix (here, the susceptibility matrix)
!          Note that wfprod of the first index must behave like a density,
!          so that it is used as generated by fourwf, while wfprod of the
!          second index will be implicitely used to make a scalar product
!          with a potential change, meaning that its complex conjugate must be
!          used. This explains the following signs...
           do ipw1=ipw2,npwdiel
             susmat(1,ipw1,jsp,ipw2,jsp)=susmat(1,ipw1,jsp,ipw2,jsp)+&
&             weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
             susmat(2,ipw1,jsp,ipw2,jsp)=susmat(2,ipw1,jsp,ipw2,jsp)+&
&             weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
           end do
         end do
       end if

       if( occopt>=3 .and. testocc==1 .and. abs(wght1)>tol12) then
!        Accumulate product of band densities by their doccde, for the
!        computation of the effect of change of Fermi level.
         do ipw=1,npwdiel
           drhode(1,ipw,jsp)=drhode(1,ipw,jsp)+wfprod(1,ipw)*wght1
           drhode(2,ipw,jsp)=drhode(2,ipw,jsp)+wfprod(2,ipw)*wght1
         end do
!        Also accumulate weighted sum of doccde
         sumdocc=sumdocc+wght1
       end if

!      End condition of metallic occupancies or extrapolation
     end if

!    End loop on spinorial components
   end do
!  End loop on iband
 end do

 call timab(752,2,tsec)
 call timab(753,1,tsec)

 ABI_FREE(cwavef)

!Stuff for parallelism (bands-FFT)
 if(mpi_enreg%paral_kgb==1) then
   call xmpi_sum(wfrspa,mpi_enreg%comm_bandfft,ier)
   if(occopt>=3) then
     call xmpi_sum(drhode(:,:,isp1:isp2),mpi_enreg%comm_bandfft,ier)
   end if
   if(extrap==1) then
     call xmpi_sum(rhoextrap,mpi_enreg%comm_bandfft,ier)
   end if
   if(occopt>=3) then
     call xmpi_sum(sumdocc,mpi_enreg%comm_bandfft,ier)
   end if
   ABI_MALLOC(susmat_mpi,(2,npwdiel,npwdiel))
   do i1=isp1,isp2
     susmat_mpi(:,:,:)=susmat(:,:,i1,:,i1)
     call xmpi_sum(susmat_mpi,mpi_enreg%comm_bandfft,ier)
     susmat(:,:,i1,:,i1)=susmat_mpi(:,:,:)/real(mpi_enreg%nproc_fft*mpi_enreg%nproc_band,dp)
   end do
   ABI_FREE(susmat_mpi)
 end if
 call timab(753,2,tsec)

!-- Wavefunctions have been generated in real space ------------------------
!-- Now, compute product of wavefunctions for different bands --------------
 call timab(754,1,tsec)
!if (occopt<3) then
 tolocc=1.0d-3
!else
!tolocc=1.0d-8
!end if
 iproc=-1

 if(nband_k>1)then
   do ibd1=1,nband_k-1
     do ibd2=ibd1+1,nband_k
       iproc=iproc+1
       if(modulo(iproc,mpi_enreg%nproc_fft*mpi_enreg%nproc_band) /= me_bandfft) cycle
!      If the occupation numbers are sufficiently different, or
!      if extrapolation is used and the corresponding factor is not zero,
!      then there is a contribution
       occdiff=occ(ibd1+bdtot_index)-occ(ibd2+bdtot_index)
       if( abs(occdiff)>tolocc      .or. &
&       ( extrap==1 .and.            &
&       ( abs(occ_deavg(ibd1)) + abs(occ_deavg(ibd2)) ) >tolocc ) &
&       ) then

         eigdiff=eigen(ibd1+bdtot_index)-eigen(ibd2+bdtot_index)
!        DEBUG
!        write(std_out,*)' susk : contribution from bands',ibd1,ibd2
!        write(std_out,*)'   occ diff =',occdiff
!        write(std_out,*)'   eig diff =',eigdiff
!        ENDDEBUG

!        Loop on spinorial components
         do ispinor=1,nspinor
           jsp=isp+ispinor-1;if (nspden_eff==1) jsp=isp

!          Store the contribution in wfraug
           do i3=1,ndiel3
             do i2=1,ndiel2
               do i1=1,ndiel1
                 wfraug(1,i1,i2,i3)=wfrspa(1,i1,i2,i3,ispinor,ibd1)*wfrspa(1,i1,i2,i3,ispinor,ibd2)&
&                 +wfrspa(2,i1,i2,i3,ispinor,ibd1)*wfrspa(2,i1,i2,i3,ispinor,ibd2)
                 wfraug(2,i1,i2,i3)=wfrspa(2,i1,i2,i3,ispinor,ibd1)*wfrspa(1,i1,i2,i3,ispinor,ibd2)&
&                 -wfrspa(1,i1,i2,i3,ispinor,ibd1)*wfrspa(2,i1,i2,i3,ispinor,ibd2)
               end do
             end do
           end do

!          Performs the Fourier Transform of the product, and store it in wfprod
           tim_fourwf=19
           call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&           1,kg_diel,kg_diel,mgfftdiel,mpi_enreg_diel,1,ngfftdiel,1,npwdiel,&
&           ndiel4,ndiel5,ndiel6,3,tim_fourwf,weight,weight)

!          In case of PAW, add compensation charge contribution
           if (usepaw==1.and.neglect_pawhat==0) then
             call pawsushat(atindx,cprj_k,gbound_diel,gylmg_diel,ibd1,ibd2,ispinor,ispinor,1,kg_diel,&
&             lmax_diel,mgfftdiel,natom,nband_k,ndiel4,ndiel5,ndiel6,&
&             ngfftdiel,npwdiel,nspinor,ntypat,0,&
&             pawang,pawtab,ph3d_diel,typat,wfprod,dummy,&
&             mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
           end if

!          Perform now the summation
           wght1=zero ; wght2=zero
           if(abs(occdiff)>tolocc) wght1= occdiff/eigdiff * two*wtk(ikpt)/ucvol
           if(extrap==1) wght2=(occ_deavg(ibd1)+occ_deavg(ibd2)) * two*wtk(ikpt)/ucvol
           weight=wght1+wght2

!          DEBUG
!          write(std_out,*)' weight =',weight
!          norm=zero
!          do ipw=1,npwdiel
!          norm=norm+wfprod(1,ipw)**2+wfprod(2,ipw)**2
!          end do
!          write(std_out,*)' norm in reciprocal space  =',norm
!          ENDDEBUG

           if (abs(weight)>tol12) then
             do ipw2=1,npwdiel
!              Only fills lower half of the matrix (here, the susceptibility matrix)
!              Note that wfprod of the first index must behave like a density,
!              so that it is used as generated by fourwf, while wfprod of the
!              second index will be implicitely used to make a scalar product
!              with a potential change, meaning that its complex conjugate must be
!              used. This explains the following signs...
               do ipw1=ipw2,npwdiel
                 susmat(1,ipw1,jsp,ipw2,jsp)=susmat(1,ipw1,jsp,ipw2,jsp)+&
&                 weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
                 susmat(2,ipw1,jsp,ipw2,jsp)=susmat(2,ipw1,jsp,ipw2,jsp)+&
&                 weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
               end do
             end do
           end if

!          End loop on spinorial components
         end do
!        End condition of different occupation numbers or extrapolation
       end if
!      End internal loop over bands
     end do
!    End external loop over bands
   end do
!  End condition of having more than one band
 end if

 call timab(754,2,tsec)
 call timab(755,1,tsec)

 if(mpi_enreg%paral_kgb==1) then
   ABI_MALLOC(susmat_mpi,(2,npwdiel,npwdiel))
   do i1=isp1,isp2
     susmat_mpi(:,:,:)=susmat(:,:,i1,:,i1)
     call xmpi_sum(susmat_mpi,mpi_enreg%comm_bandfft,ier)
     susmat(:,:,i1,:,i1)=susmat_mpi(:,:,:)
   end do
   ABI_FREE(susmat_mpi)
   ABI_FREE(band_loc)
   ABI_FREE(treat_band)
   ABI_FREE(cg_local)
   ABI_FREE(kg_k_gather_all)
 end if

 call destroy_mpi_enreg(mpi_enreg_diel)
 ABI_FREE(dummy)
 ABI_FREE(rhoaug)
 ABI_FREE(wfprod)
 ABI_FREE(wfraug)
 ABI_FREE(wfrspa)

 call timab(755,2,tsec)
 call timab(750,2,tsec)

end subroutine susk
!!***

!!****f* m_suscep_stat/suskmm
!! NAME
!! suskmm
!!
!! FUNCTION
!! Compute the contribution of one k point to the susceptibility matrix
!! from input wavefunctions, band occupations, and k point wts.
!! Include the usual sum-over-state terms, but also the
!! corrections due to the change of the Fermi level in the metallic
!! case, as well as implicit sum over higher lying conduction
!! states, thanks to the closure relation (referred to as an extrapolation).
!!
!! This routine is similar to susk, but use blocking on wavefunctions
!! to decrease memory requirements, at the expense of CPU time.
!!
!! NOTES
!! There is still room for optimization !!
!!
!! INPUTS
!!  atindx(natom)=index table for atoms
!!  bdtot_index=index for the number of the band
!!  cg(2,mcg)=wf in G space
!!  cprj_k(natom,nspinor*nband_k)= wave functions projected with non-local projectors:
!!                                 cprj_k=<p_i|Cnk> where p_i is a non-local projector.
!!  doccde(mband*nkpt*nsppol)=derivative of occupancies wrt
!!           the energy for each band and k point
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  extrap: if==1, the closure relation (an extrapolation) must be used
!!  gbound(2*mgfftdiel+8,2)=G sphere boundary for going from WF sphere to
!!      medium size FFT grid
!!  gbound_diel(2*mgfftdiel+8,2)=G sphere boundary for going from medium size
!!      FFT grid to small sphere.
!!  gylmg_diel(npwdiel,lmax_diel**2,ntypat*usepaw)= -PAW only- Fourier transform of g_l(r).Y_ml(r) shape functions
!!                                                   for dielectric matrix
!!  icg=index for cg
!!  ikpt=number of the k point
!!  isp=number of the current spin
!!  istwf_k=input option parameter that describes the storage of wfs
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  kg_k(3,npw)=coordinates of planewaves in basis sphere.
!!  lmax_diel=1+max. value of l angular momentum used for dielectric matrix
!!  mband=maximum number of bands
!!  mcg=dimension of cg
!!  mgfftdiel=maximum size of 1D FFTs, for the computation of
!!     the dielectric matrix
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell
!!  nband_k=number of bands at this k point for that spin polarization
!!  ndiel4,ndiel5,ndiel6= FFT dimensions, modified to avoid cache trashing
!!  neglect_pawhat=1 if PAW contribution from hat density (compensation charge)
!!                 has to be neglected (to be used when only an estimation of
!!                 suscep. matrix has to be evaluated, i.e. for SCF precondictioning)
!!  ngfftdiel(18)=contain all needed information about 3D FFT, for dielectric matrix,
!!    see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkpt=number of k points
!!  npwdiel=third and fifth dimension of the susmat array.
!!  npw_k=number of plane waves at this k point
!!  nspden=number of spin-density components
!!  nspden_eff=number of spin-density components actually computed in sussceptibility
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  occ(mband*nkpt*nsppol)=
!!          occupation numbers for each band (usually 2.0) at each k point
!!  occopt=option for occupancies
!!  occ_deavg(mband)=factor for extrapolation (occup. divided by an energy gap)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ph3d_diel(2,npwdiel,natom*usepaw)=3-dim structure factors, for each atom and plane wave, for dielectric matrix
!!  typat(natom)=type (integer) for each atom
!!  ucvol=unit cell volume (Bohr**3)
!!  usepaw=flag for PAW
!!  wtk(nkpt)=k point weights (they sum to 1.0)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  drhode(2,npwdiel,nspden_eff)=weighted density, needed to compute the
!!   effect of change of fermi energy
!!  rhoextrap(ndiel4,ndiel5,ndiel6,nspinor)=density-like array, needed for the
!!   extrapolation procedure.
!!  sumdocc=sum of weighted occupation numbers, needed to compute the
!!   effect of change of fermi energy
!!  susmat(2,npwdiel,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!
!! PARENTS
!!      m_suscep_stat
!!
!! CHILDREN
!!      destroy_mpi_enreg,fourwf,init_distribfft_seq,initmpi_seq,pawsushat
!!      timab
!!
!! SOURCE

subroutine suskmm(atindx,bdtot_index,cg,cprj_k,doccde,drhode,eigen,extrap,gbound,&
&  gbound_diel,gylmg_diel,icg,ikpt,isp,istwf_k,kg_diel,kg_k,&
&  lmax_diel,mband,mcg,mgfftdiel,mpi_enreg,&
&  natom,nband_k,ndiel4,ndiel5,ndiel6,neglect_pawhat,ngfftdiel,nkpt,&
&  npwdiel,npw_k,nspden,nspden_eff,nspinor,nsppol,ntypat,occ,occopt,occ_deavg,&
&  pawang,pawtab,ph3d_diel,rhoextrap,sumdocc,&
&  susmat,typat,ucvol,usepaw,wtk)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bdtot_index,extrap,icg,ikpt,isp,istwf_k,lmax_diel,mband,mcg
 integer,intent(in) :: mgfftdiel,natom,nband_k,ndiel4,ndiel5,ndiel6,neglect_pawhat
 integer,intent(in) :: nkpt,npw_k,npwdiel,nspden,nspden_eff,nspinor
 integer,intent(in) :: nsppol,ntypat,occopt,usepaw
 real(dp),intent(in) :: ucvol
 real(dp),intent(inout) :: sumdocc
 type(MPI_type),intent(in) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: atindx(natom),gbound(2*mgfftdiel+8,2)
 integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
 integer,intent(in) :: kg_diel(3,npwdiel),kg_k(3,npw_k),ngfftdiel(18)
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: cg(2,mcg),doccde(mband*nkpt*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
 real(dp),intent(in) :: gylmg_diel(npwdiel,lmax_diel**2,ntypat*usepaw)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),occ_deavg(mband)
 real(dp),intent(in) :: ph3d_diel(2,npwdiel,natom*usepaw),wtk(nkpt)
 real(dp),intent(inout) :: drhode(2,npwdiel,nspden_eff)
 real(dp),intent(inout) :: rhoextrap(ndiel4,ndiel5,ndiel6,nspinor)
 real(dp),intent(inout) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 type(pawcprj_type) :: cprj_k(natom,nspinor*nband_k*usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: comm_fft,i1,i2,i3,iband,iband_shift,iband_shift2,ibd1,ibd2,ibdshft1,ibdshft2
 integer :: iblk1,iblk2,ipw,ipw1,ipw2,ispinor,iwf,jsp,mblk
 integer :: nblk,nbnd_current,nbnd_in_blk,nbnd_in_blk1,ndiel1,ndiel2,ndiel3
 integer :: testocc,tim_fourwf
 real(dp) :: eigdiff,occdiff,tolocc,weight,wght1,wght2
 character(len=500) :: message
 type(MPI_type) :: mpi_enreg_diel
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: cwavef(:,:),dummy(:,:),rhoaug(:,:,:),wfprod(:,:)
 real(dp),allocatable :: wfraug(:,:,:,:),wfrspa1(:,:,:,:,:,:)
 real(dp),allocatable :: wfrspa2(:,:,:,:,:,:)

! *************************************************************************

 call timab(760,1,tsec)
 call timab(761,1,tsec)

!Allocations, initializations
 ndiel1=ngfftdiel(1) ; ndiel2=ngfftdiel(2) ; ndiel3=ngfftdiel(3)
 testocc=1
 ABI_MALLOC(rhoaug,(ndiel4,ndiel5,ndiel6))
 ABI_MALLOC(wfraug,(2,ndiel4,ndiel5,ndiel6))
 ABI_MALLOC(wfprod,(2,npwdiel))
 ABI_MALLOC(dummy,(2,1))

!The dielectric stuff is performed in sequential mode.
!Set mpi_enreg_diel accordingly
 call initmpi_seq(mpi_enreg_diel)
 call init_distribfft_seq(mpi_enreg_diel%distribfft,'c',ngfftdiel(2),ngfftdiel(3),'all')

 comm_fft=mpi_enreg%comm_fft

!Prepare the blocking : compute the number of blocks,
!the number of bands in each normal block,
!and the number in the first one, usually smaller.

!Consider that if the number of bands is large, there are at most 8 blocks
 nbnd_in_blk=0
 if(nband_k>=48)then
   mblk=8
   nbnd_in_blk=(nband_k-1)/mblk+1
!  If the number of bands is medium, place 6 bands per block
 else if(nband_k>=12)then
   nbnd_in_blk=6
!  Otherwise, must have at least 2 blocks
 else if(nband_k>=2)then
   mblk=2
   nbnd_in_blk=(nband_k-1)/mblk+1
 else
   write(message, '(a,a,a,i2,a,a,a)')&
&   '  The number of bands must be larger or equal to 2, in suskmm.',ch10,&
&   '  It is equal to ',nband_k,'.',ch10,&
&   '  Action : choose another preconditioner.'
   ABI_ERROR(message)
 end if

!Compute the effective number of blocks, and the number of bands in
!the first block.
 nblk=(nband_k-1)/nbnd_in_blk+1
 nbnd_in_blk1=nband_k-(nblk-1)*nbnd_in_blk

!DEBUG
!write(std_out,*)' suskmm : nband_k,nblk,nbnd_in_blk,nbnd_in_blk1 '
!write(std_out,*)nband_k,nblk,nbnd_in_blk,nbnd_in_blk1
!stop
!ENDDEBUG

!wfrspa1 will contain the wavefunctions of the slow sampling (iblk1)
 ABI_MALLOC(wfrspa1,(2,ndiel4,ndiel5,ndiel6,nspinor,nbnd_in_blk))
!wfrspa2 will contain the wavefunctions of the rapid sampling (iblk2)
 ABI_MALLOC(wfrspa2,(2,ndiel4,ndiel5,ndiel6,nspinor,nbnd_in_blk))

 ABI_MALLOC(cwavef,(2,npw_k))

 call timab(761,2,tsec)

!First loop over blocks
 do iblk1=1,nblk

   call timab(762,1,tsec)

!  Initialisation
   if(iblk1==1)then

     nbnd_current=nbnd_in_blk1
     iband_shift=0
!    Loop over bands to fft and store Fourier transform of wavefunction
     do iband=1,nbnd_current
!      Loop on spinorial components
       do ispinor=1,nspinor
         iwf=(ispinor-1)*npw_k+(iband-1)*npw_k*nspinor+icg
!        Obtain Fourier transform in fft box
         tim_fourwf=21
         cwavef(:,1:npw_k)=cg(:,1+iwf:npw_k+iwf)
         call fourwf(1,rhoaug,cwavef,dummy,wfraug,gbound,gbound,&
&         istwf_k,kg_k,kg_k,mgfftdiel,mpi_enreg_diel,1,ngfftdiel,npw_k,1,ndiel4,ndiel5,ndiel6,&
&         0,tim_fourwf,weight,weight)
         wfrspa1(:,:,:,:,ispinor,iband)=wfraug(:,:,:,:)
       end do
     end do

   else

!    The Fourier transform of wavefunctions have already been obtained
     nbnd_current=nbnd_in_blk
     iband_shift=nbnd_in_blk1+(iblk1-2)*nbnd_in_blk

   end if

!  Loop over bands of this block, to generate band-diagonal
   do iband=1,nbnd_current

!    Loop on spinorial components
     do ispinor=1,nspinor
       jsp=isp+ispinor-1;if (nspden_eff==1) jsp=isp

       if( (occopt>=3 .and. testocc==1) .or. extrap==1 )then
!        In the case of metallic occupation, or if the extrapolation
!        over higher bands is included, must compute the
!        Fourier transform of the density of each band, then
!        generate the part of the susceptibility matrix due
!        varying occupation numbers.
         weight=-two*occ_deavg(iband+iband_shift)*wtk(ikpt)/ucvol
         do i3=1,ndiel3
           do i2=1,ndiel2
             do i1=1,ndiel1
               wfraug(1,i1,i2,i3)=wfrspa1(1,i1,i2,i3,ispinor,iband)**2&
&               +wfrspa1(2,i1,i2,i3,ispinor,iband)**2
               wfraug(2,i1,i2,i3)=zero
             end do
           end do
!          If extrapolation, accumulate density in real space
           if(extrap==1.and.usepaw==0)then
             do i2=1,ndiel2
               do i1=1,ndiel1
                 rhoextrap(i1,i2,i3,ispinor)=rhoextrap(i1,i2,i3,ispinor)+weight*wfraug(1,i1,i2,i3)
               end do
             end do
           end if
         end do

!        In case of PAW, add compensation charge contribution
         if (usepaw==1.and.extrap==1.and.neglect_pawhat==0) then
           call pawsushat(atindx,cprj_k,gbound_diel,gylmg_diel,iband,iband,ispinor,ispinor,1,kg_diel,&
&           lmax_diel,mgfftdiel,natom,nband_k,ndiel4,ndiel5,ndiel6,&
&           ngfftdiel,npwdiel,nspinor,ntypat,1,&
&           pawang,pawtab,ph3d_diel,typat,dummy,wfraug,&
&           mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
           rhoextrap(:,:,:,ispinor)=rhoextrap(:,:,:,ispinor)+weight*wfraug(1,:,:,:)
         end if

!        Performs the Fourier Transform of the density of the band,
!        and store it in wfprod
         tim_fourwf=31
         call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&         1,kg_diel,kg_diel,&
&         mgfftdiel,mpi_enreg_diel,1,ngfftdiel,1,npwdiel,ndiel4,ndiel5,ndiel6,3,tim_fourwf,weight,weight)
!        In case of PAW, add compensation charge contribution if not already done
         if (usepaw==1.and.extrap==0.and.neglect_pawhat==0) then
           call pawsushat(atindx,cprj_k,gbound_diel,gylmg_diel,iband,iband,1,1,1,kg_diel,&
&           lmax_diel,mgfftdiel,natom,nband_k,ndiel4,ndiel5,ndiel6,&
&           ngfftdiel,npwdiel,nspinor,ntypat,0,&
&           pawang,pawtab,ph3d_diel,typat,wfprod,dummy,&
&           mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
         end if

!        Perform now the summation of terms related to direct change of eigenvalues
!        or extrapolation over higher bands
         wght1=zero ; wght2=zero
         if(occopt>=3 .and. testocc==1) wght1=doccde(iband+iband_shift+bdtot_index)*wtk(ikpt)/ucvol
         if(extrap==1) wght2=two*occ_deavg(iband+iband_shift)*wtk(ikpt)/ucvol
         weight=wght1+wght2

         if (abs(weight)>tol12) then
           do ipw2=1,npwdiel
!            Only fills lower half of the matrix (here, the susceptibility matrix)
!            Note that wfprod of the first index must behave like a density,
!            so that it is used as generated by fourwf, while wfprod of the
!            second index will be implicitely used to make a scalar product
!            with a potential change, meaning that its complex conjugate must be
!            used. This explains the following signs...
             do ipw1=ipw2,npwdiel
               susmat(1,ipw1,jsp,ipw2,jsp)=susmat(1,ipw1,jsp,ipw2,jsp)+&
&               weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
               susmat(2,ipw1,jsp,ipw2,jsp)=susmat(2,ipw1,jsp,ipw2,jsp)+&
&               weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
             end do
           end do
         end if

         if( occopt>=3 .and. testocc==1 .and. abs(wght1)>tol12) then
!          Accumulate product of band densities by their doccde, for the
!          computation of the effect of change of Fermi level.
           do ipw=1,npwdiel
             drhode(1,ipw,jsp)=drhode(1,ipw,jsp)+wfprod(1,ipw)*wght1
             drhode(2,ipw,jsp)=drhode(2,ipw,jsp)+wfprod(2,ipw)*wght1
           end do
!          Also accumulate weighted sum of doccde
           sumdocc=sumdocc+wght1
         end if

!        End condition of metallic occupancies or extrapolation
       end if

!      End loop on spinorial components
     end do
!    End loop on iband
   end do

   call timab(762,2,tsec)

!  -- Compute now off-band-diagonal terms ------------------------------------
!  -- Compute product of wavefunctions for different bands, inside the block -

   call timab(763,1,tsec)

!  if (occopt<3) then
   tolocc=1.0d-3
!  else
!  tolocc=1.0d-8
!  end if

   if(nbnd_current>1)then
     do ibd1=1,nbnd_current-1
       ibdshft1=ibd1+iband_shift
       do ibd2=ibd1+1,nbnd_current
         ibdshft2=ibd2+iband_shift

!        If the occupation numbers are sufficiently different, or
!        if extrapolation is used and the corresponding factor is not zero,
!        then there is a contribution
         occdiff=occ(ibdshft1+bdtot_index)-occ(ibdshft2+bdtot_index)
         if( abs(occdiff)>tolocc      .or. &
&         ( extrap==1 .and.            &
&         ( abs(occ_deavg(ibdshft1)) + abs(occ_deavg(ibdshft2)) ) >tolocc ) &
&         ) then

           eigdiff=eigen(ibdshft1+bdtot_index) - eigen(ibdshft2+bdtot_index)

!          Loop on spinorial components
           do ispinor=1,nspinor
             jsp=isp+ispinor-1;if (nspden_eff==1) jsp=isp

!            Store the contribution in wfraug
             do i3=1,ndiel3
               do i2=1,ndiel2
                 do i1=1,ndiel1
                   wfraug(1,i1,i2,i3)=wfrspa1(1,i1,i2,i3,ispinor,ibd1)*wfrspa1(1,i1,i2,i3,ispinor,ibd2)&
&                   +wfrspa1(2,i1,i2,i3,ispinor,ibd1)*wfrspa1(2,i1,i2,i3,ispinor,ibd2)
                   wfraug(2,i1,i2,i3)=wfrspa1(2,i1,i2,i3,ispinor,ibd1)*wfrspa1(1,i1,i2,i3,ispinor,ibd2)&
&                   -wfrspa1(1,i1,i2,i3,ispinor,ibd1)*wfrspa1(2,i1,i2,i3,ispinor,ibd2)
                 end do
               end do
             end do

!            Performs the Fourier Transform of the product, and store it in wfprod
             tim_fourwf=32
             call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&             1,kg_diel,kg_diel, mgfftdiel,mpi_enreg_diel,1,ngfftdiel,1,npwdiel,&
&             ndiel4,ndiel5,ndiel6,3,tim_fourwf,weight,weight)

!            In case of PAW, add compensation charge contribution
             if (usepaw==1.and.neglect_pawhat==0) then
               call pawsushat(atindx,cprj_k,gbound_diel,gylmg_diel,ibd1,ibd2,ispinor,ispinor,1,kg_diel,&
&               lmax_diel,mgfftdiel,natom,nband_k,ndiel4,ndiel5,ndiel6,&
&               ngfftdiel,npwdiel,nspinor,ntypat,0,&
&               pawang,pawtab,ph3d_diel,typat,wfprod,dummy,&
&               mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
             end if

!            Perform now the summation
             wght1=zero ; wght2=zero
             if(abs(occdiff)>tolocc) wght1= occdiff/eigdiff * two*wtk(ikpt)/ucvol
             if(extrap==1) wght2=(occ_deavg(ibdshft1)+occ_deavg(ibdshft2)) * two*wtk(ikpt)/ucvol
             weight=wght1+wght2

             if (abs(weight)>tol12) then
               do ipw2=1,npwdiel
!                Only fills lower half of the matrix (here, the susceptibility matrix)
!                Note that wfprod of the first index must behave like a density,
!                so that it is used as generated by fourwf, while wfprod of the
!                second index will be implicitely used to make a scalar product
!                with a potential change, meaning that its complex conjugate must be
!                used. This explains the following signs...
                 do ipw1=ipw2,npwdiel
                   susmat(1,ipw1,jsp,ipw2,jsp)=susmat(1,ipw1,jsp,ipw2,jsp)+&
&                   weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
                   susmat(2,ipw1,jsp,ipw2,jsp)=susmat(2,ipw1,jsp,ipw2,jsp)+&
&                   weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
                 end do
               end do
             end if

!            End loop on spinorial components
           end do
!          End condition of different occupation numbers or extrapolation
         end if
!        End internal loop over bands
       end do
!      End external loop over bands
     end do
!    End condition of having more than one band
   end if

!  Loop on secondary block, with fast varying index, in decreasing order.
   if(iblk1/=nblk)then
     do iblk2=nblk,iblk1+1,-1
       iband_shift2=nbnd_in_blk1+(iblk2-2)*nbnd_in_blk

!      Loop over bands to fft and store Fourier transform of wavefunction
       iband_shift2=nbnd_in_blk1+(iblk2-2)*nbnd_in_blk
       do iband=1,nbnd_in_blk
!        Loop on spinorial components
         do ispinor=1,nspinor
           iwf=(ispinor-1)*npw_k+(iband+iband_shift2-1)*npw_k*nspinor+icg

!          Obtain Fourier transform in fft box
           tim_fourwf=22
           cwavef(:,1:npw_k)=cg(:,1+iwf:npw_k+iwf)
           call fourwf(1,rhoaug,cwavef,dummy,wfraug,gbound,gbound,&
&           istwf_k,kg_k,kg_k,mgfftdiel,mpi_enreg_diel,1,ngfftdiel,npw_k,1,&
&           ndiel4,ndiel5,ndiel6,0,tim_fourwf,weight,weight)
           wfrspa2(:,:,:,:,ispinor,iband)=wfraug(:,:,:,:)
         end do
       end do

       do ibd1=1,nbnd_current
         ibdshft1=ibd1+iband_shift
         do ibd2=1,nbnd_in_blk
           ibdshft2=ibd2+iband_shift2

!          If the occupation numbers are sufficiently different, or
!          if extrapolation is used and the corresponding factor is not zero,
!          then there is a contribution
           occdiff=occ(ibdshft1+bdtot_index)-occ(ibdshft2+bdtot_index)
           if( abs(occdiff)>tolocc      .or. &
&           ( extrap==1 .and.            &
&           ( abs(occ_deavg(ibdshft1)) + abs(occ_deavg(ibdshft2)) ) >tolocc ) &
&           ) then

             eigdiff=eigen(ibdshft1+bdtot_index) - eigen(ibdshft2+bdtot_index)

!            Loop on spinorial components
             do ispinor=1,nspinor
               jsp=isp+ispinor-1;if (nspden_eff==1) jsp=isp

!              Store the contribution in wfraug
               do i3=1,ndiel3
                 do i2=1,ndiel2
                   do i1=1,ndiel1
                     wfraug(1,i1,i2,i3)=wfrspa1(1,i1,i2,i3,ispinor,ibd1)*wfrspa2(1,i1,i2,i3,ispinor,ibd2)&
&                     +wfrspa1(2,i1,i2,i3,ispinor,ibd1)*wfrspa2(2,i1,i2,i3,ispinor,ibd2)
                     wfraug(2,i1,i2,i3)=wfrspa1(2,i1,i2,i3,ispinor,ibd1)*wfrspa2(1,i1,i2,i3,ispinor,ibd2)&
&                     -wfrspa1(1,i1,i2,i3,ispinor,ibd1)*wfrspa2(2,i1,i2,i3,ispinor,ibd2)
                   end do
                 end do
               end do

!              Performs the Fourier Transform of the product, and store it in wfprod
               tim_fourwf=32
               call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&               1,kg_diel,kg_diel,mgfftdiel,mpi_enreg_diel,1,ngfftdiel,1,npwdiel,&
&               ndiel4,ndiel5,ndiel6,3,tim_fourwf,weight,weight)

!              In case of PAW, add compensation charge contribution
               if (usepaw==1.and.neglect_pawhat==0) then
                 call pawsushat(atindx,cprj_k,gbound_diel,gylmg_diel,ibd1,ibdshft2,ispinor,ispinor,1,kg_diel,&
&                 lmax_diel,mgfftdiel,natom,nband_k,ndiel4,ndiel5,ndiel6,&
&                 ngfftdiel,npwdiel,nspinor,ntypat,0,&
&                 pawang,pawtab,ph3d_diel,typat,wfprod,dummy,&
&                 mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
               end if

!              Perform now the summation
               wght1=zero ; wght2=zero
               if(abs(occdiff)>tolocc) wght1= occdiff/eigdiff * two*wtk(ikpt)/ucvol
               if(extrap==1) wght2=(occ_deavg(ibdshft1)+occ_deavg(ibdshft2)) * two*wtk(ikpt)/ucvol
               weight=wght1+wght2

               if (abs(weight)>tol12) then
                 do ipw2=1,npwdiel
!                  Only fills lower half of the matrix (here, the susceptibility matrix)
!                  Note that wfprod of the first index must behave like a density,
!                  so that it is used as generated by fourwf, while wfprod of the
!                  second index will be implicitely used to make a scalar product
!                  with a potential change, meaning that its complex conjugate must be
!                  used. This explains the following signs...
                   do ipw1=ipw2,npwdiel
                     susmat(1,ipw1,jsp,ipw2,jsp)=susmat(1,ipw1,jsp,ipw2,jsp)+&
&                     weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
                     susmat(2,ipw1,jsp,ipw2,jsp)=susmat(2,ipw1,jsp,ipw2,jsp)+&
&                     weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
                   end do
                 end do
               end if

!              End loop on spinorial components
             end do
!            End condition of different occupation numbers or extrapolation
           end if
!          End internal loop over bands
         end do
!        End external loop over bands
       end do
!      End loop on bloks
     end do

!    Finish the loop on blok with iblk2=iblk1+1, so can use the
!    FFTd wavefunctions for the next iblk1.
     do iband=1,nbnd_in_blk
       wfrspa1(:,:,:,:,1:nspinor,iband)=wfrspa2(:,:,:,:,1:nspinor,iband)
     end do

!    End condition of iblk1/=nblk
   end if

   call timab(763,2,tsec)

!  End loop on iblk1
 end do

!DEBUG
!write(std_out,*)' suskmm : exit '
!do ipw1=1,npwdiel
!write(std_out,*)ipw1,susmat(1,ipw1,1,ipw1,1),susmat(2,ipw1,1,ipw1,1)
!end do
!write(std_out,*)' suskmm : end of susmat '
!stop
!ENDDEBUG

 call destroy_mpi_enreg(mpi_enreg_diel)
 ABI_FREE(cwavef)
 ABI_FREE(dummy)
 ABI_FREE(rhoaug)
 ABI_FREE(wfprod)
 ABI_FREE(wfraug)
 ABI_FREE(wfrspa1)
 ABI_FREE(wfrspa2)

 call timab(760,2,tsec)

end subroutine suskmm
!!***

end module m_suscep_stat
!!***
