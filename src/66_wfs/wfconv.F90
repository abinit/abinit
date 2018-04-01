!{\src2tex{textfont=tt}}
!!****f* ABINIT/wfconv
!! NAME
!! wfconv
!!
!! FUNCTION
!! This subroutine treats the wavefunctions for one k point,
!! and converts them to other parameters.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2018 ABINIT group (XG,TD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ceksp2=if 1, center the output sphere of pw on Gamma; if 0, on each k-point (usual).
!!  cg1(2,mcg1)=wavefunction array
!!  debug= if 1, print some messages ; otherwise, 0.
!!  ecut1=kinetic energy cutoffs for basis sphere 1 (hartree)
!!  ecut2=kinetic energy cutoff beyond which the coefficients of wf2 vanish (Ha)
!!  ecut2_eff=kinetic energy cut-off for basis sphere 2 (hartree)
!!  eig_k1(mband1*(2*mband1)**formeig)=eigenvalues
!!  exchn2n3d=if 1, n2 and n3 are exchanged
!!  formeig option (format of the eigenvalues and eigenvector) :
!!   0 => ground-state format (initialisation of
!!        eigenvectors with random numbers, vector of eigenvalues)
!!   1 => respfn format (initialisation of
!!        eigenvectors with 0 s, hermitian matrix of eigenvalues)
!!  gmet1(3,3)=reciprocal space metric (bohr^-2) for input wf
!!  gmet2(3,3)=reciprocal space metric (bohr^-2) for output wf
!!  icg1=shift to be given to the location of the data in the array cg1
!!  icg2=shift to be given to the location of the data in the array cg2
!!  ikpt1=number of the k point actually treated (input wf numbering)
!!  ikpt10=number of the k point previously treated (input wf numbering)
!!  ikpt2=number of the k point actually treated (output numbering)
!!  indkk(nkpt2*sppoldbl,6)=describe k point number of kptns1 that allows to
!!   generate wavefunctions closest to given kpt2 (and possibly isppol2=2)
!!   indkk(:,1)=k point number of kpt1
!!   indkk(:,2)=symmetry operation to be applied to kpt1, to give kpt1a
!!    (if 0, means no symmetry operation, equivalent to identity )
!!   indkk(:,3:5)=shift in reciprocal space to be given to kpt1a,
!!    to give kpt1b, that is the closest to kpt2.
!!   indkk(:,6)=1 if time-reversal was used to generate kpt1a from kpt1, 0 otherwise
!!  inplace= if 0, cg1 and cg2 are different in the calling routine,
!!           if 1, cg1 and cg2 are identical (they have the same memory location)
!!    This is also true for the pairs (eig_k1,eig_k2) and (occ_k1,occ_k2)
!!  isppol2=spin variable for output wavefunctions
!!  istwfk1(nkpt1)=input parameter that describes the storage of wfs in set1
!!  istwfk2(nkpt2)=input parameter that describes the storage of wfs in set2
!!  kg1(3,mpw1)=dimensionless coords of G vecs in basis sphere at k point (input wf)
!!  kg2(3,mpw2)=dimensionless coords of G vecs in basis sphere at k point (output wf)
!!  kptns1(3,nkpt1)=k point set for input wavefunctions
!!  kptns2(3,nkpt2)=k point set for output wavefunctions
!!  mband1=dimension of eig_k1 and occ_k1 arrays
!!  mband2=dimension of eig_k2 and occ_k2 arrays
!!  mcg1=dimension of cg1 array (at least npw1*nspinor1*nbd1)
!!  mcg2=dimension of cg2 array (at least npw2*nspinor2*nbd2)
!!  mpi_enreg1=informations about MPI parallelization for set 1
!!  mpi_enreg2=informations about MPI parallelization for set 2
!!  mpw1=dimension of kg1, can be set to 0 if not needed
!!  mpw2=dimension of kg2, can be set to 0 if not needed
!!  nbd1=number of bands contained in cg1,eig_k1,occ_k1 at this k-point - spin (at input)
!!  nbd2=number of bands contained in cg2,eig_k2,occ_k2 at this k-point - spin (at output)
!!  ngfft1(18)=all needed information about 3D FFT, for input wavefunctions
!!  ngfft2(18)=all needed information about 3D FFT, for output wavefunctions
!!             see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkpt1=number of k points for input wavefunctions
!!  nkpt2=number of k points for output wavefunctions
!!  npw1=number of planewaves for input wavefunctions
!!  npw2=number of planewaves for output wavefunctions
!!  nspinor1=number of spinors for input wavefunctions
!!  nspinor2=number of spinors for output wavefunctions
!!  nsym=number of symmetry elements in space group
!!  occ_k1(mband1)=occupation numbers
!!  optorth=1 if the WFs are orthogonalized before leaving the routine
!!  randalg=1 if "good" (but non-portable) random numbers should be used, 0 for compatibility
!!  restart=if 2, conversion between wavefunctions
!!          if 1, direct restart is allowed (see hdr_check.f)
!!  rprimd2(3,3)=dimensional primitive translations for real space (bohr)
!!   needed only for the spinor rotation
!!  sppoldbl= if 1, no doubling of the number if spins thanks to antiferromagn
!!    if 2, deduce nsppol=2 from nsppol=1, using Shubnikov symmetries
!!  symrel(3,3,nsym)=symmetry operations in real space in terms
!!   of primitive translations
!!  tnons(3,nsym)=nonsymmorphic translations for symmetry operations
!!
!! OUTPUT
!!  cg2(2,mcg2)=wavefunction array
!!  eig_k2(mband2*(2*mband2)**formeig)=eigenvalues
!!  occ_k2(mband2)=occupation (completed with zeros)
!!
!! SIDE EFFECTS
!! Input/Output:
!!  ikpt10=at input, number of the k point previously treated (input wf numbering)
!!     (if this is the first call for the present k point set, ikpt10 should be 0)
!!         at output, number of the k point just treated (input wf numbering)
!!  kg1, kg2, npw1 and npw2 should not be modified by kpgsph (TD).
!!
!! NOTES
!! Note that this routine can make an in-place conversion
!! (see the input variable "inplace"),
!! if cg1 and cg2 are equal, as well as the pairs (icg1,icg2),
!! (eig_k1,eig_k2),(occ_k1,occ_k2) and (mband1,mband2)
!!
!! It can also be used to fill or to initialize wavefunctions
!! at one k point
!! (filling with random numbers or 0''s, according to the value
!! of formeig), if the input number of bands (nbd1) is 0.
!! In the latter case, one should use the same values of input
!! wavefunction parameters
!! than for output wavefunction parameters, except nbd1.
!!
!! The input parameters are indexed with 1, the output parameters
!! are indexed with 2.
!!
!! Some of the arguments are arrays dimensioned with nkpt1 or nkpt2.
!! Note that for these, only the elements for ikpt1 or ikpt2 will be used.
!!
!! The number of input bands must already be minimal at the input.
!! This means, when input and output nspinor are equal : nbd1<nbd2
!! When the two nspinor differ, one must have nbd1/nspinor1<nbd2/nspinor2
!!
!! PARENTS
!!      newkpt,wfsinp
!!
!! CHILDREN
!!      cg_envlop,getph,getspinrot,kpgsph,mati3inv,ph1d3d,pw_orthon,sphere
!!      sphereboundary,timab,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wfconv(ceksp2,cg1,cg2,debug,ecut1,ecut2,ecut2_eff,&
& eig_k1,eig_k2,exchn2n3d,formeig,gmet1,gmet2,icg1,icg2,&
& ikpt1,ikpt10,ikpt2,indkk,inplace,isppol2,istwfk1,istwfk2,&
& kg1,kg2,kptns1,kptns2,mband1,mband2,mcg1,mcg2,mpi_enreg1,mpi_enreg2,&
& mpw1,mpw2,nbd1,nbd2,ngfft1,ngfft2,nkpt1,nkpt2,npw1,npw2,nspinor1,nspinor2,&
& nsym,occ_k1,occ_k2,optorth,randalg,restart,rprimd2,sppoldbl,symrel,tnons)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

 use m_time,     only : timab
 use m_fftcore,  only : kpgsph, sphere, sphereboundary
 use m_cgtools,  only : cg_envlop
 use m_kg,       only : ph1d3d, getph

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfconv'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_66_wfs, except_this_one => wfconv
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ceksp2,debug,exchn2n3d,formeig,icg1,icg2,ikpt1
 integer,intent(in) :: ikpt2,inplace,isppol2,mband1,mband2,mcg1,mcg2,mpw1,mpw2
 integer,intent(in) :: nbd1,nbd2,nkpt1,nkpt2,nspinor1,nspinor2,nsym
 integer,intent(in) :: optorth,randalg,restart,sppoldbl
 integer,intent(inout) :: ikpt10,npw1,npw2
 real(dp),intent(in) :: ecut1,ecut2,ecut2_eff
 type(MPI_type),intent(inout) :: mpi_enreg1,mpi_enreg2
!arrays
 integer,intent(in) :: indkk(nkpt2*sppoldbl,6),istwfk1(nkpt1),istwfk2(nkpt2)
 integer,intent(in) :: ngfft1(18),ngfft2(18),symrel(3,3,nsym)
 integer,intent(inout) :: kg1(3,mpw1),kg2(3,mpw2)
 real(dp),intent(in) :: gmet1(3,3),gmet2(3,3),kptns1(3,nkpt1),kptns2(3,nkpt2)
 real(dp),intent(in) :: rprimd2(3,3),tnons(3,nsym)
 real(dp),intent(inout) :: cg1(2,mcg1),cg2(2,mcg2)
 real(dp),intent(inout) :: eig_k1(mband1*(2*mband1)**formeig)
 real(dp),intent(inout) :: eig_k2(mband2*(2*mband2)**formeig),occ_k1(mband1)
 real(dp),intent(inout) :: occ_k2(mband2)

!Local variables ------------------------------
!scalars
 integer,parameter :: nkpt_max=50,tobox=1,tosph=-1
 integer :: conv_tnons,convert,fftalg,fold1,fold2,foldim,foldre,i1,i2,iband
 integer :: iband_first,iband_last,icgmod,ierr,index,ipw
 integer :: ispinor,ispinor1,ispinor2,ispinor_first,ispinor_last
 integer :: istwf10_k,istwf1_k,istwf2_k,isym,itimrev,jsign
 integer :: mgfft1,mgfft2,n1,n2,n3,n4,n5,n6
 integer :: nbremn,npwtot,nspinor_index,nspinor1_this_proc,nspinor2_this_proc
 integer :: order,ortalgo,seed
 real(dp) :: ai,ar,arg,bi,br,eig_tmp,spinrots,spinrotx,spinroty,spinrotz
 character(len=500) :: message
 integer, parameter :: int64 = selected_int_kind(18)
 !arrays
 integer :: atindx(1),identity(3,3),ngfft_now(18),no_shift(3),shiftg(3)
 integer :: symm(3,3),symrel_conv(3,3)
 integer,allocatable :: gbound1(:,:),gbound2(:,:)
 real(dp) :: kpoint1(3),kpoint2_sph(3),phktnons(2,1),spinrot(4),tnons_conv(3),tsec(2)
 real(dp),allocatable :: cfft(:,:,:,:),dum(:,:),phase1d(:,:),phase3d(:,:)
 real(dp),allocatable :: wavef1(:,:),wavef2(:,:),wavefspinor(:,:)

! *************************************************************************

 mgfft1=maxval(ngfft1(1:3))
 mgfft2=maxval(ngfft2(1:3))
 if(.false.)write(std_out,*)occ_k1 ! just to keep occ_k1 as an argument before resolving the issue of its transfer

 if(nspinor1/=1 .and. nspinor1/=2)then
   write(message,'(a,i0)')'The argument nspinor1 must be 1 or 2, while it is nspinor1 = ',nspinor1
   MSG_BUG(message)
 end if

 if(nspinor2/=1 .and. nspinor2/=2)then
   write(message,'(a,i0)')' The argument nspinor2 must be 1 or 2, while it is nspinor2=',nspinor2
   MSG_BUG(message)
 end if

 if(nspinor1==2 .and. mod(nbd1,2)/=0)then
   write(message,'(a,i0)')' When nspinor1 is 2, nbd1 must be even, while it is nbd1 = ',nbd1
   MSG_BUG(message)
 end if

 if(nspinor2==2 .and. mod(nbd2,2)/=0)then
   write(message,'(a,i0)')'  When nspinor2 is 2, nbd2 must be even, while it is nbd2=',nbd2
   MSG_BUG(message)
 end if

 if(nbd1/nspinor1>nbd2/nspinor2)then
   write(message, '(3a,2i6,3a,2i6,a)' )&
&   'In wfconv, the nbd/nspinor ratio cannot decrease. However,',ch10,&
&   'the initial quantities are nbd1,nspinor1=',nbd1,nspinor1,', and',ch10,&
&   'the requested final quantities are nbd2,nspinor2=',nbd2,nspinor2,'.'
   MSG_BUG(message)
 end if

 ngfft_now(1:3)=ngfft1(1:3)
 ngfft_now(8:18)=ngfft1(8:18)
!This line is the reason why ngfft_now has to be introduced
 ngfft_now(7)=101
 ngfft_now(4:6)=ngfft_now(1:3)
 n1=ngfft_now(1) ; n2=ngfft_now(2) ; n3=ngfft_now(3)
 n4=ngfft_now(4) ; n5=ngfft_now(5) ; n6=ngfft_now(6)
 fftalg=ngfft_now(7)

!Parallelization over spinors management
 nspinor1_this_proc=max(1,nspinor1/mpi_enreg1%nproc_spinor)
 nspinor2_this_proc=max(1,nspinor2/mpi_enreg2%nproc_spinor)

!In order to generate IN PLACE new wfs from old wfs, the loop
!over bands and spinors must be done in one direction or the other,
!depending on npw1 and npw2, nspinor1 and nspinor2.
!If nspinor1=1 and nspinor2=2 , note that one will generate
!from nbd1 states of npw1 coefficients,
!2*nbd1 states of 2*npw2 coefficients. nbd1 cancels in comparing
!these expressions, but nspinor2 appears squared.
!The same line of thought works for the case nspinor1=2 and nspinor2=1
 order=1
 iband_first=1   ; iband_last=nbd1
 ispinor_first=1 ; ispinor_last=nspinor1
 if(nspinor1==2 .and. nspinor2==1)then
   order=2 ; iband_last=nbd1-1 ; ispinor_last=1
 end if
!Here, reverse the order if needed
 if( npw2*nspinor2**2 > npw1*nspinor1**2 )then
   order=-order
   iband_first=iband_last       ; iband_last=1
   ispinor_first=ispinor_last   ; ispinor_last=1
 end if

 kpoint1(:)=kptns1(:,ikpt1)
 istwf1_k=istwfk1(ikpt1)

 kpoint2_sph(:)=0.0_dp
 if(ceksp2==0)kpoint2_sph(:)=kptns2(:,ikpt2)
 istwf2_k=istwfk2(ikpt2)

!DEBUG
!write(std_out,*)'ecut1,ecut2_eff=',ecut1,ecut2_eff
!write(std_out,*)'gmet1,gmet2=',gmet1,gmet2
!write(std_out,*)'kpoint1,kpoint2_sph=',kpoint1,kpoint2_sph
!write(std_out,*)'nspinor1,nspinor2',nspinor1,nspinor2
!write(std_out,*)'istwf1_k,istwf2_k=',istwf1_k,istwf2_k
!write(std_out,*)'nbd1,tol8=',nbd1,tol8
!ENDDEBUG

!Determine whether it will be needed to convert the existing
!wavefunctions, or simply to complete them.

 convert=0
 if(nbd1/=0)then
   if(abs(ecut2_eff-ecut1)>tol8)convert=convert+1
   if(sum(abs(gmet2(:,:)-gmet1(:,:)))>tol8)convert=convert+2
   if(sum(abs(kpoint2_sph(:)-kpoint1(:)))>tol8)convert=convert+4
   if(nspinor2/=nspinor1)convert=convert+8
   if(istwf2_k/=istwf1_k)convert=convert+16
 end if

!This is a supplementary check
 if(restart==1 .and. convert/=0)then
   MSG_BUG('Restart==1 and convert/=0 are exclusive')
 end if

!Determine whether symmetries must be used
 conv_tnons=0
 no_shift(:)=0
 identity(:,:)=0
 identity(1,1)=1 ; identity(2,2)=1 ; identity(3,3)=1
 isym=indkk(ikpt2+(sppoldbl-1)*(isppol2-1)*nkpt2,2)
!DEBUG
!write(std_out,*)' wfconv : isym=',isym
!ENDDEBUG
 itimrev=indkk(ikpt2+(sppoldbl-1)*(isppol2-1)*nkpt2,6)
 if(isym/=0)then
   symrel_conv(:,:)=symrel(:,:,isym)
   call mati3inv(symrel_conv,symm)
   shiftg(:)=indkk(ikpt2+(sppoldbl-1)*(isppol2-1)*nkpt2,3:5)
   tnons_conv(:)=tnons(:,isym)
   if(sum(tnons_conv(:)**2)>tol8)then
!    Need to compute phase factors associated with nonsymmorphic translations.
     conv_tnons=1
     ABI_ALLOCATE(phase3d,(2,npw1))
     ABI_ALLOCATE(phase1d,(2,(2*n1+1)+(2*n2+1)+(2*n3+1)))
!    Although the routine getph is originally written for
!    atomic phase factors, it does precisely what we want
     atindx(1)=1
     call getph(atindx,1,n1,n2,n3,phase1d,tnons_conv)
   end if
   if(nspinor1==2 .and. nspinor2==2)then
!    Compute rotation in spinor space
     call getspinrot(rprimd2,spinrot,symrel_conv)
   end if
 else
   shiftg(:)=0
   symm(:,:)=identity(:,:)
   spinrot(:)=zero
   spinrot(1)=one
 end if
 if(itimrev/=0)then
   symm(:,:)=-symm(:,:)
 end if

!DEBUG
!write(std_out,'(a,i3,2x,3i3,2x,9i3)')' wfconv : isym,shiftg,symm=',isym,shiftg,symm
!write(std_out,*)' wfconv : ecut2_eff,ecut1=',ecut2_eff,ecut1
!write(std_out,*)' wfconv : istwf1_k,istwf2_k=',istwf1_k,istwf2_k
!write(std_out,*)' wfconv : kpoint1(:),kpoint2_sph(:)=',&
!& kpoint1(:),kpoint2_sph(:)
!write(std_out,*)' wfconv : nspinor1,nspinor2=',nspinor1,nspinor2
!ENDDEBUG

!if (mpi_enreg1%fft_option_lob==0) mpi_enreg1%fft_option_lob=1
!if (mpi_enreg2%fft_option_lob==0) mpi_enreg2%fft_option_lob=1

 if (restart==2.and.(convert/=0.or.(nbd2/nspinor2>nbd1/nspinor1.and.formeig==0))) then
!  kg2 is needed both for FFT grid conversion and for envlop
!  Choose the center of the sphere : either gamma, or each k-point
   kpoint2_sph(:)=0.0_dp
   if(ceksp2==0)kpoint2_sph(:)=kptns2(:,ikpt2)
   istwf2_k=istwfk2(ikpt2)
   call kpgsph(ecut2_eff,exchn2n3d,gmet2,0,ikpt2,istwf2_k,kg2,kpoint2_sph,1,mpi_enreg2,mpw2,npw2)
 end if

 if(convert/=0)then
   istwf10_k=0
   if(ikpt10/=0)istwf10_k=istwfk1(ikpt10)

!  Only need G sphere if different from last time
   if ( ikpt1/=ikpt10 .or. istwf1_k/=istwf10_k ) then

     call kpgsph (ecut1,exchn2n3d,gmet1,0,ikpt1,istwf1_k,kg1,kpoint1,1,mpi_enreg1,mpw1,npw1)
     if (debug>0) then
       write(message, '(a,f8.3,a,a,3f8.5,a,a,i3,a,3(a,3es16.8,a),a,3i4,a,i5,a)' )&
&       ' wfconv: called kpgsph with ecut1=',ecut1,ch10,&
&       '  kpt1=',kptns1(1:3,ikpt1),ch10,&
&       '  istwf1_k=',istwf1_k,ch10,&
&       '  gmet1= ',gmet1(1:3,1),ch10,&
&       '         ',gmet1(1:3,2),ch10,&
&       '         ',gmet1(1:3,3),ch10,&
&       '  ngfft=',ngfft_now(1:3),' giving npw1=',npw1,'.'
       call wrtout(std_out,message,'PERS')
     end if
     ikpt10 = ikpt1
     istwf10_k=istwf1_k
   end if

   if(conv_tnons==1)then
     arg=two_pi*(kpoint1(1)*tnons_conv(1)+ kpoint1(2)*tnons_conv(2)+ kpoint1(3)*tnons_conv(3) )
     phktnons(1,1)=cos(arg)
     phktnons(2,1)=sin(arg)
!    Convert 1D phase factors to 3D phase factors exp(i 2 pi (k+G).tnons )
     call ph1d3d(1,1,kg1,1,1,npw1,n1,n2,n3,phktnons,phase1d,phase3d)
   end if

   ABI_ALLOCATE(cfft,(2,n4,n5,n6))
   ABI_ALLOCATE(wavef1,(2,npw1))
   ABI_ALLOCATE(wavef2,(2,npw2))
   if(nspinor1==2 .and. nspinor2==2) then
     ABI_ALLOCATE(wavefspinor,(2,2*npw2))
   end if
   ABI_ALLOCATE(gbound1,(2*mgfft1+8,2))
   ABI_ALLOCATE(gbound2,(2*mgfft2+8,2))
   call sphereboundary(gbound1,istwf1_k,kg1,mgfft1,npw1)
   call sphereboundary(gbound2,istwf2_k,kg2,mgfft2,npw2)

!  Take old wf from sphere->box, the new from box->sphere
!  One pays attention not to have a problem of erasing data when replacing
!  a small set of coefficient by a large set, or the reverse.
!  This is the reason of the use of order, _first and _last variables,
!  defined earlier.
   nspinor_index=mpi_enreg1%me_spinor+1
   do iband=iband_first,iband_last,order
     do ispinor1=ispinor_first,ispinor_last,order
       ispinor=ispinor1
       if (mpi_enreg1%paral_spinor==1) then
         if (ispinor1==nspinor_index) then
           ispinor=1
         else
           if (nspinor1==2.and.nspinor2==2) wavefspinor(:,(ispinor1-1)*npw2+1:ispinor1*npw2)=zero
           cycle
         end if
       end if

!      Copy input wf
       i1=(ispinor-1)*npw1+(iband-1)*nspinor1_this_proc*npw1+icg1
       wavef1(:,1:npw1)=cg1(:,i1+1:i1+npw1)

!      Make symmetry-induced conversion, if needed (translation part)
       if(conv_tnons==1)then
!$OMP PARALLEL DO PRIVATE(ai,ar)
         do ipw=1,npw1
           ar=phase3d(1,ipw)*wavef1(1,ipw)-phase3d(2,ipw)*wavef1(2,ipw)
           ai=phase3d(2,ipw)*wavef1(1,ipw)+phase3d(1,ipw)*wavef1(2,ipw)
           wavef1(1,ipw)=ar
           wavef1(2,ipw)=ai
         end do
       end if

!      Take into account time-reversal symmetry, if needed, in the scalar case
       if(itimrev==1 .and. (nspinor1==1 .or. nspinor2==1))then
!$OMP PARALLEL DO 
         do ipw=1,npw1
           wavef1(2,ipw)=-wavef1(2,ipw)
         end do
       end if

!      DEBUG
!      write(std_out,*)' wfconv : before sphere, isym,ispinor=',isym,ispinor
!      write(std_out,*)' no_shift,identity=',no_shift,identity
!      write(std_out,*)' shiftg,symm=',shiftg,symm
!      stop
!      This debugging sequence is an attempt to rotate spinors,
!      and works indeed for test13, when symmetry 9 is used ...
!      if(isym==9 .and. ispinor==1)then
!      write(std_out,*)' wfconv : gives a 120 degree rotation to first component'
!      do ipw=1,npw1
!      ar=-            half*wavef1(1,ipw)-sqrt(three)*half*wavef1(2,ipw)
!      ai= sqrt(three)*half*wavef1(1,ipw)-            half*wavef1(2,ipw)
!      wavef1(1,ipw)=ar
!      wavef1(2,ipw)=ai
!      end do
!      end if
!      ENDDEBUG

!      Convert wf, and also include the symmetry operation and shiftg.
       call sphere(wavef1,1,npw1,cfft,n1,n2,n3,n4,n5,n6,kg1,istwf1_k,tobox,&
&       mpi_enreg1%me_g0,no_shift,identity,one)

       call sphere(wavef2,1,npw2,cfft,n1,n2,n3,n4,n5,n6,kg2,istwf2_k,tosph,&
&       mpi_enreg2%me_g0,shiftg,symm,one)

       if(nspinor2==1 )then
         i2=(ispinor-1)*npw2+(iband-1)*nspinor2_this_proc*npw2+icg2
         cg2(:,i2+1:i2+npw2)=wavef2(:,1:npw2)
       else if(nspinor1==2.and.nspinor2==2)then
!        Will treat this case outside of the ispinor loop
         i2=(ispinor1-1)*npw2
         wavefspinor(:,i2+1:i2+npw2)=wavef2(:,1:npw2)
       else if(nspinor1==1 .and. nspinor2==2)then
!        The number of bands is doubled, and the number of coefficients
!        is doubled also
         if (mpi_enreg2%paral_spinor==0) then
           i2=(iband-1)*nspinor2_this_proc*nspinor2_this_proc*npw2+icg2
           cg2(:,i2+1:i2+npw2)=wavef2(:,1:npw2)
           cg2(:,i2+npw2+1:i2+2*npw2)=zero
           cg2(:,i2+2*npw2+1:i2+3*npw2)=zero
           cg2(:,i2+3*npw2+1:i2+4*npw2)=wavef2(:,1:npw2)
         else
           i2=(iband-1)*nspinor2_this_proc*npw2+icg2
           if (nspinor_index==1) then
             cg2(:,i2+1:i2+npw2)=wavef2(:,1:npw2)
             cg2(:,i2+npw2+1:i2+2*npw2)=zero
           else
             cg2(:,i2+1:i2+npw2)=zero
             cg2(:,i2+npw2+1:i2+2*npw2)=wavef2(:,1:npw2)
           end if
         end if
       end if
     end do ! ispinor=ispinor_first,ispinor_last,order

     if(nspinor1==2.and.nspinor2==2)then
!      Take care of possible parallelization over spinors
       if (mpi_enreg2%paral_spinor==1) then
         call xmpi_sum(wavefspinor,mpi_enreg2%comm_spinor,ierr)
       end if
!      Take care of time-reversal symmetry, if needed
       if(itimrev==1)then
!        Exchange spin-up and spin-down
!        Make complex conjugate of one component,
!        and change sign of other component
!$OMP PARALLEL DO PRIVATE(ipw,ar,ai) SHARED(wavefspinor,npw2)
         do ipw=1,npw2
!          Here, change sign of real part
           ar=-wavefspinor(1,ipw)
           ai= wavefspinor(2,ipw)
           wavefspinor(1,ipw)= wavefspinor(1,npw2+ipw)
!          Here, change sign of imaginary part
           wavefspinor(2,ipw)=-wavefspinor(2,npw2+ipw)
           wavefspinor(1,npw2+ipw)=ar
           wavefspinor(2,npw2+ipw)=ai
         end do
       end if ! itimrev==1

!      Rotation in spinor space
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(npw2,spinrot,wavefspinor)
       spinrots=spinrot(1)
       spinrotx=spinrot(2)
       spinroty=spinrot(3)
       spinrotz=spinrot(4)
!$OMP DO
       do ipw=1,npw2
         ar=wavefspinor(1,ipw)
         ai=wavefspinor(2,ipw)
         br=wavefspinor(1,npw2+ipw)
         bi=wavefspinor(2,npw2+ipw)
         wavefspinor(1,ipw)     = spinrots*ar-spinrotz*ai +spinroty*br-spinrotx*bi
         wavefspinor(2,ipw)     = spinrots*ai+spinrotz*ar +spinroty*bi+spinrotx*br
         wavefspinor(1,npw2+ipw)=-spinroty*ar-spinrotx*ai +spinrots*br+spinrotz*bi
         wavefspinor(2,npw2+ipw)=-spinroty*ai+spinrotx*ar +spinrots*bi-spinrotz*br
       end do
!$OMP END DO
!$OMP END PARALLEL

!      Save wavefunction
       i2=(iband-1)*nspinor2_this_proc*npw2+icg2
       if (mpi_enreg2%paral_spinor==0) then
         cg2(:,i2     +1:i2+  npw2)=wavefspinor(:,1:npw2)
         cg2(:,i2+npw2+1:i2+2*npw2)=wavefspinor(:,npw2+1:2*npw2)
       else
         if (nspinor_index==1) then
           cg2(:,i2+1:i2+npw2)=wavefspinor(:,1:npw2)
         else
           cg2(:,i2+1:i2+npw2)=wavefspinor(:,npw2+1:2*npw2)
         end if
       end if
     end if ! nspinor1==2 .and. nspinor2==2

   end do

!  Take care of copying eig and occ when nspinor increases or decreases
   if(nspinor1==1.and.nspinor2==2)then
     if(formeig==0)then
!      Note the reverse order, needed in case inplace=1
       do iband=nbd1,1,-1
!        use eig_tmp to avoid bug on ifort10.1 x86_64
         eig_tmp=eig_k1(iband)
         eig_k2(2*iband-1:2*iband)=eig_tmp
!        occ_tmp=occ_k1(iband)*0.5_dp
!        occ_k2(2*iband-1:2*iband )=occ_tmp
       end do
     else
       call wrtout(std_out,' wfconv: not yet coded, formeig=1!',"COLL")
     end if
   end if
   if(nspinor1==2 .and. nspinor2==1)then
     if(formeig==0)then
       do iband=1,nbd1
!        use eig_tmp to avoid bug on ifort10.1 x86_64
         eig_tmp=eig_k1(2*iband-1)
         eig_k2(iband)=eig_tmp
!        occ_tmp=occ_k1(2*iband-1)*2.0_dp
!        occ_k2(iband)=occ_tmp
       end do
     else
       call wrtout(std_out,' wfconv: not yet coded, formeig=1!',"COLL")
     end if
   end if

   ABI_DEALLOCATE(cfft)
   ABI_DEALLOCATE(gbound1)
   ABI_DEALLOCATE(gbound2)
   ABI_DEALLOCATE(wavef1)
   ABI_DEALLOCATE(wavef2)
   if(nspinor1==2 .and. nspinor2==2) then
     ABI_DEALLOCATE(wavefspinor)
   end if

 else if(convert==0)then

   if(inplace==0)then
!    Must copy cg, eig and occ if not in-place while convert==0
!    Note that npw1=npw2, nspinor1=nspinor2
     cg2(:,1+icg2:npw1*nspinor1_this_proc*nbd1+icg2)=&
&     cg1(:,1+icg1:npw1*nspinor1_this_proc*nbd1+icg1)
     eig_k2(:)=eig_k1(:)
!    occ_k2(:)=occ_k1(:)
   end if

 end if ! End of if convert/=0

 if(conv_tnons==1) then 
   ABI_DEALLOCATE(phase1d)
   ABI_DEALLOCATE(phase3d)
 end if


!If not enough bands, complete with random numbers or zeros
 if(nbd2/nspinor2>nbd1/nspinor1)then
   if(formeig==0)then

!    Ground state wf and eig case
     eig_k2((nbd1/nspinor1)*nspinor2+1:nbd2)=huge(0.0_dp)/10.0_dp
     occ_k2((nbd1/nspinor1)*nspinor2+1:nbd2)=0.0_dp
     index=(nbd1/nspinor1)*nspinor2*npw2*nspinor2_this_proc

!    Initialisation of wavefunctions
!    One needs to initialize wfs in such a way to avoid symmetry traps,
!    and to avoid linear dependencies between wavefunctions
!    No need for a difference for different k points and/or spin-polarization

     if (mpi_enreg1%paral_kgb == 1) then
       npwtot=npw2
       call timab(539,1,tsec)
       call xmpi_sum(npwtot, mpi_enreg1%comm_bandfft, ierr)
       call timab(539,2,tsec)
     end if

     do iband=(nbd1/nspinor1)*nspinor2+1,nbd2
       do ispinor2=1,nspinor2_this_proc
         ispinor=ispinor2;if (nspinor2_this_proc/=nspinor2) ispinor=mpi_enreg2%me_spinor+1
         jsign=1;if (ispinor==2) jsign=-1
         
         do ipw=1,npw2
           index=index+1
!          Different seed for different planewave and band
!          DEBUG seq==par
!          if(.false.) then
!          ENDDEBUG seq==par

           if ( mpi_enreg2%paral_kgb /= 1) then
             seed=(iband-1)*npw2*nspinor2 + (ispinor-1)*npw2 + ipw
           else
             seed=jsign*(iband*(kg2(1,ipw)*npwtot*npwtot + kg2(2,ipw)*npwtot + kg2(3,ipw)))
           end if
           
           if(randalg == 0) then
!            For portability, use only integer numbers
!            The series of couples (fold1,fold2) is periodic with a period of
!            3x5x7x11x13x17x19x23x29x31, that is, larger than 2**32, the largest integer*4
!            fold1 is between 0 and 34, fold2 is between 0 and 114. As sums of five
!            uniform random variables, their distribution is close to a gaussian
             fold1=mod(seed,3)+mod(seed,5)+mod(seed,7)+mod(seed,11)+mod(seed,13)
             fold2=mod(seed,17)+mod(seed,19)+mod(seed,23)+mod(seed,29)+mod(seed,31)
!            The gaussian distributions are folded, in order to be back to a uniform distribution
!            foldre is between 0 and 20, foldim is between 0 and 18
             foldre=mod(fold1+fold2,21)
             foldim=mod(3*fold1+2*fold2,19)
             cg2(1,index+icg2)=dble(foldre)
             cg2(2,index+icg2)=dble(foldim)
           else
             ! (AL) Simple linear congruential generator from
             ! numerical recipes, modulo'ed and 64bit'ed to avoid
             ! overflows (NAG doesn't like overflows, even though
             ! they are perfectly legitimate here). Then, we get some
             ! lowest order bits and sum them, as the previous
             ! generator, to get quasi-normal numbers.
             ! This is clearly suboptimal and might cause problems,
             ! but at least it doesn't seem to create linear
             ! dependencies and local minima like the previous one.
             ! it's not trivial to generate good reproductible random
             ! numbers in parallel. Patches welcome !
             ! Note a fun fortran fact : MOD simply ignores 64 bits integer
             ! and casts them into 32bits, so we use MODULO.
             fold1 = modulo(1664525_int64 * seed  + 1013904223_int64, 2147483648_int64)
             fold2 = modulo(1664525_int64 * fold1 + 1013904223_int64, 2147483648_int64)
             fold1=modulo(fold1,3)+modulo(fold1,5)+modulo(fold1,7)+modulo(fold1,11)+modulo(fold1,13)
             fold2=modulo(fold2,3)+modulo(fold2,5)+modulo(fold2,7)+modulo(fold2,11)+modulo(fold2,13)
             cg2(1,index+icg2)=dble(fold1)/34-0.5
             cg2(2,index+icg2)=dble(fold2)/34-0.5
           end if
         end do
       end do
       
!      XG030513 : MPIWF need to impose cg to zero when at Gamma
!      Time-reversal symmetry for k=gamma impose zero imaginary part at G=0
!      XG : I do not know what happens for spin-orbit here :
       if(istwf2_k==2 .and. mpi_enreg2%me_g0==1) then 
         cg2(2,1+(iband-1)*npw2*nspinor2_this_proc+icg2)=zero
       end if
     end do

!    Multiply with envelope function to reduce kinetic energy
     icgmod=icg2+npw2*nspinor2_this_proc*(nbd1/nspinor1)
     nbremn=nbd2-nbd1
     call cg_envlop(cg2,ecut2,gmet2,icgmod,kg2,kpoint2_sph,mcg2,nbremn,npw2,nspinor2_this_proc)

     if(ikpt2<=nkpt_max)then
       write(message,'(3(a,i6))')' wfconv:',nbremn,' bands initialized randomly with npw=',npw2,', for ikpt=',ikpt2
       call wrtout(std_out,message,'PERS')
     end if

   else if(formeig==1)then

!    For response function, put large numbers in the remaining of the
!    eigenvalue array (part of it was already filled in calling routine)
!    WARNING : Change of nspinor not yet coded
     eig_k2(1+2*nbd1*nbd2 : 2*nbd2*nbd2)=huge(0.0_dp)/10.0_dp
!    Initialisation of wfs with 0 s
     index=npw2*nbd1*nspinor2_this_proc
     do iband=nbd1+1,nbd2
       do ipw=1,npw2*nspinor2_this_proc
         index=index+1
         cg2(:,index+icg2)=zero
       end do
     end do

     if(ikpt2<=nkpt_max)then
       nbremn=nbd2-nbd1
       write(message,'(a,i6,a,i7,a,i4)')' wfconv :',nbremn,' bands set=0 with npw=',npw2,', for ikpt=',ikpt2
       call wrtout(std_out,message,'PERS')
     end if

   end if ! End of initialisation to 0
 end if

!Orthogonalize GS wfs
 !if (.False.) then
 if (optorth==1.and.formeig==0.and.mpi_enreg2%paral_kgb/=1) then
   ABI_ALLOCATE(dum,(2,0))
   ortalgo=0 !;ortalgo=3
   call pw_orthon(icg2,0,istwf2_k,mcg2,0,npw2*nspinor2_this_proc,nbd2,ortalgo,dum,0,cg2,&
&   mpi_enreg2%me_g0,mpi_enreg2%comm_bandspinorfft)
   ABI_DEALLOCATE(dum)
 end if

end subroutine wfconv
!!***
