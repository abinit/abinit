!{\src2tex{textfont=tt}}
!!****f* ABINIT/suscep_stat
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
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (XG,AR,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!!      vtorho
!!
!! CHILDREN
!!      destroy_mpi_enreg,fftpac,init_distribfft_seq,initmpi_seq,pawcprj_alloc
!!      pawcprj_free,pawcprj_get,pawcprj_mpi_allgather,pawgylmg,ph1d3d
!!      sphereboundary,susk,suskmm,symg,symrhg,timab,xmpi_sum,zhpev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine suscep_stat(atindx,atindx1,cg,cprj,dielar,dimcprj,doccde,&
&  eigen,gbound_diel,gprimd,irrzondiel,istwfk,kg,&
&  kg_diel,lmax_diel,&
&  mband,mcg,mcprj,mgfftdiel,mkmem,mpi_enreg,mpw,natom,nband,&
&  neglect_pawhat,nfftdiel,ngfftdiel,nkpt,npwarr,&
&  npwdiel,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&
&  pawang,pawtab,phnonsdiel,ph1ddiel,rprimd,&
&  susmat,symafm,symrel,tnons,typat,ucvol,unpaw,usecprj,usepaw,usetimerev,&
&  wtk,ylmdiel)

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_profiling_abi

 use m_pawang,  only : pawang_type
 use m_pawtab,  only : pawtab_type
 use m_pawcprj, only : pawcprj_type, pawcprj_alloc, &
&                      pawcprj_get, pawcprj_mpi_allgather, pawcprj_free
 use m_mpinfo,  only : destroy_mpi_enreg
 use m_kg,      only : ph1d3d

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'suscep_stat'
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_51_manage_mpi
 use interfaces_52_fft_mpi_noabirule
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_65_paw
 use interfaces_67_common
 use interfaces_77_suscep, except_this_one => suscep_stat
!End of the abilint section

 implicit none

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
 integer :: paral_kgb_diel,spaceComm,t1,t2,testocc
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
   MSG_ERROR(message)
 end if

 if (mpi_enreg%paral_spinor==1) then
   message = ' not yet allowed for parallelization over spinors !'
   MSG_ERROR(message)
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
 ABI_ALLOCATE(occ_deavg,(mband))
 if(occopt>=3) then
   ABI_ALLOCATE(drhode,(2,npwdiel,nspden_eff))
 else
   ABI_ALLOCATE(drhode,(0,0,0))
 end if
 if(extrap==1) then
   ABI_ALLOCATE(rhoextrap,(ndiel4,ndiel5,ndiel6,nspinor))
 else
   ABI_ALLOCATE(rhoextrap,(0,0,0,0))
 end if

!zero the susceptibility matrix and other needed quantities
 susmat(:,:,:,:,:)=zero
 if(occopt>=3)then
   drhode(:,:,:)=zero
   sumdocc=zero
 end if

!PAW additional initializations
 if (usepaw==1) then
   ABI_ALLOCATE(gylmg_diel,(npwdiel,lmax_diel**2,ntypat))
   ABI_ALLOCATE(ph3d_diel,(2,npwdiel,natom))
   if (neglect_pawhat==0) then
     ABI_ALLOCATE(phkxred_diel,(2,natom))
     ABI_ALLOCATE(kpg_dum,(0,0))
     kpt_diel(1:3,1)=zero;phkxred_diel(1,:)=one;phkxred_diel(2,:)=zero;nkpg_diel=0
!    write(std_out,*) ' lmax_diel ', lmax_diel
     call pawgylmg(gprimd,gylmg_diel,kg_diel,kpg_dum,kpt_diel,lmax_diel,nkpg_diel,npwdiel,ntypat,pawtab,ylmdiel)
     call ph1d3d(1,natom,kg_diel,natom,natom,npwdiel,ndiel1,ndiel2,ndiel3,phkxred_diel,ph1ddiel,ph3d_diel)
     ABI_DEALLOCATE(phkxred_diel)
     ABI_DEALLOCATE(kpg_dum)
   else
     gylmg_diel=zero;ph3d_diel=one
   end if
 else
   ABI_ALLOCATE(gylmg_diel,(0,0,0))
   ABI_ALLOCATE(ph3d_diel,(0,0,0))
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


     ABI_ALLOCATE(gbound,(2*mgfftdiel+8,2))
     ABI_ALLOCATE(kg_k,(3,npw_k))

     if (usepaw==1) then
       ABI_DATATYPE_ALLOCATE(cprj_k,(natom,my_nspinor*nband_k))
       if (neglect_pawhat==0) then
         call pawcprj_alloc(cprj_k,0,dimcprj)
         if (mpi_enreg%nproc_band==1) then
           call pawcprj_get(atindx1,cprj_k,cprj,natom,1,ibg,ikpt,iorder_cprj,isp,&
&           mband,mkmem,natom,nband_k,nband_k,my_nspinor,nsppol,unpaw,&
&           mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         else
           nband_loc=nband_k/mpi_enreg%nproc_band
           ABI_DATATYPE_ALLOCATE(cprj_loc,(natom,my_nspinor*nband_loc))
           call pawcprj_alloc(cprj_loc,0,dimcprj)
           call pawcprj_get(atindx1,cprj_loc,cprj,natom,1,ibg,ikpt,iorder_cprj,isp,&
&           mband/mpi_enreg%nproc_band,mkmem,natom,nband_loc,nband_loc,my_nspinor,nsppol,unpaw,&
&           mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
           call pawcprj_mpi_allgather(cprj_loc,cprj_k,natom,my_nspinor*nband_loc,dimcprj,0,&
&           mpi_enreg%nproc_band,mpi_enreg%comm_band,ierr,rank_ordered=.true.)
           call pawcprj_free(cprj_loc)
           ABI_DATATYPE_DEALLOCATE(cprj_loc)
         end if
       else
         !call pawcprj_nullify(cprj_k)
       end if
     else
       ABI_DATATYPE_ALLOCATE(cprj_k,(0,0))
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
&       npwdiel,npw_k,nspden,nspden_eff,nspinor,nsppol,ntypat,occ,occopt,occ_deavg,paral_kgb_diel,&
&       pawang,pawtab,ph3d_diel,rhoextrap,sumdocc,&
&       susmat,typat,ucvol,usepaw,wtk)
     end if

     call timab(743,2,tsec)

     ABI_DEALLOCATE(gbound)
     ABI_DEALLOCATE(kg_k)

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
     ABI_DATATYPE_DEALLOCATE(cprj_k)

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
     ABI_ALLOCATE(rhoextrr,(nfftdiel,nspden_tmp))
     ABI_ALLOCATE(rhoextrg,(2,nfftdiel))
     if (nspden==1.and.nspinor==2) rhoextrap(:,:,:,1)=rhoextrap(:,:,:,1)+rhoextrap(:,:,:,2)

     do ispinor=1,min(nspinor,nspden)
       jsp=isp+ispinor-1

       call fftpac(1,mpi_enreg_diel,1,ndiel1,ndiel2,ndiel3,ndiel4,ndiel5,ndiel6,&
&       ngfftdiel,rhoextrr(:,1),rhoextrap(:,:,:,ispinor),1)

!      Generate the density in reciprocal space, and symmetrize it
!      (note symrhg also make the reverse FFT, to get symmetrized density;
!      this is useless here, and should be made an option)
       call symrhg(1,gprimd,irrzondiel,mpi_enreg_diel,nfftdiel,nfftdiel,ngfftdiel,&
&       nspden_tmp,1,nsym,paral_kgb_diel,phnonsdiel,rhoextrg,rhoextrr,rprimd,symafm,symrel)

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
     ABI_DEALLOCATE(rhoextrg)
     ABI_DEALLOCATE(rhoextrr)

     call timab(744,2,tsec)

   end if

!  End loop over spins ---------------------------------------------------------
 end do

 ABI_DEALLOCATE(occ_deavg)
 ABI_DEALLOCATE(rhoextrap)
 ABI_DEALLOCATE(gylmg_diel)
 ABI_DEALLOCATE(ph3d_diel)
 !end if

 call destroy_mpi_enreg(mpi_enreg_diel)

!-- Stuff for parallelism --------------------------------------------------------
!---------------------------------------------------------------------------------

 if(xmpi_paral==1)then
   call timab(746,1,tsec)
   ABI_ALLOCATE(sussum,(2*npwdiel*nspden*npwdiel*nspden))
!  Recreate full susmat on all proc.
!  This should be coded more efficiently,
!  since half of the matrix is still empty, and
!  it is spin-diagonal.
   sussum(:)=reshape(susmat(:,:,:,:,:),(/2*npwdiel*nspden*npwdiel*nspden/))
   call xmpi_sum(sussum,spaceComm,ierr)
   susmat(:,:,:,:,:)=reshape(sussum(:),(/2,npwdiel,nspden,npwdiel,nspden/))
   ABI_DEALLOCATE(sussum)
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
 ABI_ALLOCATE(tmrev_g,(npwdiel))
 ABI_ALLOCATE(sym_g,(npwdiel,nsym))
 ABI_ALLOCATE(phdiel,(2,npwdiel,nsym))
 call symg(kg_diel,npwdiel,nsym,phdiel,sym_g,symrel,tmrev_g,tnons)

!Impose spatial symmetries to the spin-diagonal susceptibility matrix
 ABI_ALLOCATE(suswk,(2,npwdiel,npwdiel,nspden_eff))
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
 ABI_DEALLOCATE(suswk)


!--  Add contribibution to susceptibility due to change of Fermi level -----------
!---------------------------------------------------------------------------------

 if (occopt>=3.and.testocc==1) then

!  Impose spatial symmetries to drhode
   ABI_ALLOCATE(drhode_wk,(2,npwdiel,nspden_eff))
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
   ABI_DEALLOCATE(drhode_wk)

 end if
 !if (occopt>=3)  then
 ABI_DEALLOCATE(drhode)
 !end if


!--- Impose the time-reversal symmetry to the susceptibility matrix --------------
!---------------------------------------------------------------------------------

 if (usetimerev==1) then
   ABI_ALLOCATE(suswk,(2,npwdiel,npwdiel,1))

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
   ABI_DEALLOCATE(suswk)
 end if

 ABI_DEALLOCATE(phdiel)
 ABI_DEALLOCATE(sym_g)
 ABI_DEALLOCATE(tmrev_g)


!-- The full susceptibility matrix is computed -----------------------------------
!-- Now, eventually diagonalize it and stop --------------------------------------
!---------------------------------------------------------------------------------

!Must turn on this flag to make the diagonalisation
 diag=0
 if(diag==1)then

   npwsp=npwdiel*nspden_eff
   ABI_ALLOCATE(sush,(npwsp*(npwsp+1)))
   ABI_ALLOCATE(susvec,(2,npwsp,npwsp))
   ABI_ALLOCATE(eig_diel,(npwsp))
   ABI_ALLOCATE(zhpev1,(2,2*npwsp-1))
   ABI_ALLOCATE(zhpev2,(3*npwsp-2))
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

   ABI_DEALLOCATE(sush)
   ABI_DEALLOCATE(susvec)
   ABI_DEALLOCATE(eig_diel)
   ABI_DEALLOCATE(zhpev1)
   ABI_DEALLOCATE(zhpev2)
   MSG_ERROR("Stopping here!")
 end if

 call timab(747,2,tsec)
 call timab(740,2,tsec)

end subroutine suscep_stat
!!***
