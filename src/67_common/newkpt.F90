!{\src2tex{textfont=tt}}
!!****f* ABINIT/newkpt
!! NAME
!! newkpt
!!
!! FUNCTION
!! This subroutine writes a starting guess for wave function (set 2)
!! It performs a "zero order" interpolation, ie simply
!! searches the nearest available k-point.
!! The data (set 1) associated with this point is either
!! read from a disk file (with a random access reading routine),
!! or input as argument.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR, ZL, AR, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  ceksp2=if 1, center the sphere of pw on Gamma; if 0, on each k-point.
!!  doorth=1 to do orthogonalization
!!  debug=>0 for debugging output
!!  ecut1=kinetic energy cutoffs for basis sphere 1 (hartree)
!!  ecut2=kinetic energy cutoffs beyond which the coefficients of wf2 vanish (Ha)
!!  ecut2_eff=kinetic energy cut-off for basis sphere 2 (hartree)
!!  exchn2n3d=if 1, n2 and n3 are exchanged
!!  fill=if 1, fill the supplementary bands ; if 0, reduce the number of bands
!!             Note : must have fill/=0 in the parallel execution
!!  formeig=if 0, GS format for wfs, eig and occ ; if 1, RF format.
!!  gmet1(3,3), gmet2(3,3)=reciprocal space metrics (bohr^-2)
!!  headform1=header format (might be needed to read the block of wfs)
!!  indkk(nkpt2*sppoldbl,6)=describe k point number of kptns1 that allows to
!!   generate wavefunctions closest to given kpt2 (and possibly isppol2=2)
!!   indkk(:,1)=k point number of kpt1
!!   indkk(:,2)=symmetry operation to be applied to kpt1, to give kpt1a
!!    (if 0, means no symmetry operation, equivalent to identity )
!!   indkk(:,3:5)=shift in reciprocal space to be given to kpt1a,
!!    to give kpt1b, that is the closest to ikpt2.
!!   indkk(:,6)=1 if time-reversal was used to generate kpt1a from kpt1, 0 otherwise
!!  iout=unit number for output file
!!  ireadwf=if 0, no reading of disk wavefunction file (random or 0.0 initialisation)
!!  istwfk1(nkpt1)=input parameter that describes the storage of wfs in set1
!!  istwfk2(nkpt2)=input parameter that describes the storage of wfs in set2
!!  kg2(3,mpw2*mkmem2)=dimensionless coords of G vecs in basis sphere at k point
!!  kptns1(3,nkpt1), kptns2(3,nkpt2)=k point sets (reduced coordinates)
!!  mband2= maximum number of bands of the output wavefunctions
!!  mcg=dimension of the cg array
!!   In case mkmem2/=0, all the output data must find their place in cg,
!!    so that mcg must be at least Sum(ikpt,isppol) [npw*nspinor*nband](ikpt,isppol)
!!    where these data are related to the output parameters
!!   In case mkmem1/=0, the same is true, for the input parameters,
!!    however, the maximum number of bands that will be read
!!    will be at most (mband2/nspinor2)*nspinor1
!!   In case mkmem1==0 and mkmem2==0, one must have at least mpw*nspinor*mband
!!    for BOTH the input and output parameters, taking into account the
!!    maximal number of band to be read, described above.
!!   In case mkmem1/=0 and mkmem2/=0, it is expected that the input cg array
!!    is organised using the output parameters nkpt2, nband2 ...
!!    This is needed, in order to use the same pointer.
!!  mkmem1= if 0, the input wf, eig, occ are available from disk
!!  mkmem2= if 0, the output wf, eig, occ must be written onto disk
!!  mpi_enreg1=informations about MPI parallelization, for the input wf file
!!  mpi_enreg2=informations about MPI parallelization, for the output wf file
!!  mpw1=maximum allowed number of planewaves at any k, for the input wf file
!!  mpw2=maximum allowed number of planewaves at any k, for the output wf file
!!  my_nkpt2= number of k points for the output wf file, handled by current processus
!!  nband1(nkpt1*nsppol1)=number of bands, at each k point, on disk
!!  nband2(nkpt2*nsppol2)=desired number of bands at each k point
!!  ngfft1(18)=all needed information about 3D FFT, for the input wf file
!!  ngfft2(18)=all needed information about 3D FFT, for the output wf file
!!             see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkpt1, nkpt2=number of k points in each set
!!  npwarr1(nkpt1)=array holding npw for each k point (input wf file).
!!  npwarr2(nkpt2)=array holding npw for each k point (output wf file).
!!  nspinor1,nspinor2=number of spinorial components of the wavefunctions
!!   for each wf file (input or output)
!!  nsppol1=1 for unpolarized, 2 for spin-polarized, input wf file
!!  nsppol2=1 for unpolarized, 2 for spin-polarized, output wf file
!!  nsym=number of symmetry elements in space group
!!  optorth= 1 if the WFS have to be orthogonalized; 0 otherwise
!!  prtvol=control print volume and debugging
!!  randalg=1 if "good" (but non-portable) random numbers should be used, 0 for compatibility
!!  restart= if 2, conversion between wavefunctions
!!           if 1, direct restart is allowed (see hdr_check.f)
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  sppoldbl= if 1, no doubling of the number if spins thanks to antiferromagn
!!    if 2, deduce nsppol=2 from nsppol=1, using Shubnikov symmetries
!!  symrel(3,3,nsym)=symmetry operations in real space in terms
!!   of primitive translations
!!  tnons(3,nsym)=nonsymmorphic translations for symmetry operations
!!  unkg2=unit number for storage of basis sphere data: stores indirect
!!   indexing array and integer coordinates for all planewaves in basis
!!   sphere for each k point being considered (kptns2 set)
!!  wffinp=structure info of input wf file unit number
!!  wffout=structure info of output wf file unit number
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!     The following arrays are input if mkmem1/=0, otherwise their input
!!     values are taken from disk, and are output if mkmem2/=0, otherwise
!!     their output values are written on disk.
!!     The location of the block for a given spin-k point at input MUST
!!     be the same as the location of the corresponding spin-k point at output.
!!  cg(2,mcg)=complex wf array
!!  eigen(mband2*(2*mband2)**formeig *nkpt2*nsppol2)=
!!    eigenvalues (input or init to large number for GS or init to 0.0 for RF), (Ha)
!!  occ(mband2*nkpt2*nsppol2)=occupation (input or init to 0.0)  NOT USED NOW
!!
!! NOTES
!! * When reading from disk, it is expected that the next record of
!! the wffinp%unwff disk unit is the first record of the first wavefunction block.
!!
!! * When the data is input as argument, it is assumed that the
!! data for each spin- k wavefunction block is located at the proper
!! corresponding location of the output array (this is to be described).
!!
!! * The information is pumped onto an fft box for the conversion.
!! This allows for changing the number of plane waves.
!!
!! * In the present status of this routine, occ is not output.
!!
!! PARENTS
!!      inwffil
!!
!! CHILDREN
!!      pareigocc,prmat,randac,rdnpw,rwwf,timab,wfconv,wffreadskipk,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine newkpt(ceksp2,cg,debug,ecut1,ecut2,ecut2_eff,eigen,exchn2n3d,fill,&
&                  formeig,gmet1,gmet2,headform1,indkk,iout,ireadwf,&
&                  istwfk1,istwfk2,kg2,kptns1,kptns2,mband2,mcg,mkmem1,mkmem2,&
&                  mpi_enreg1,mpi_enreg2,mpw1,mpw2,my_nkpt2,nband1,nband2,&
&                  ngfft1,ngfft2,nkpt1,nkpt2,npwarr1,npwarr2,nspinor1,nspinor2,&
&                  nsppol1,nsppol2,nsym,occ,optorth,prtvol,randalg,restart,rprimd,&
&                  sppoldbl,symrel,tnons,unkg2,wffinp,wffout)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_wffile
 use m_xmpi

 use m_pptools,    only : prmat
 use m_occ,        only : pareigocc

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'newkpt'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_56_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_66_wfs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ceksp2,debug,exchn2n3d,fill,formeig,headform1,iout
 integer,intent(in) :: ireadwf,mband2,mcg,mkmem1,mkmem2,mpw1,mpw2,my_nkpt2,nkpt1,nkpt2
 integer,intent(in) :: nspinor1,nspinor2,nsppol1,nsppol2,nsym,optorth,prtvol,restart
 integer,intent(in) :: randalg,sppoldbl,unkg2
 real(dp),intent(in) :: ecut1,ecut2,ecut2_eff
 type(MPI_type),intent(inout) :: mpi_enreg1,mpi_enreg2
 type(wffile_type),intent(inout) :: wffinp,wffout
!arrays
 integer,intent(in) :: indkk(nkpt2*sppoldbl,6),istwfk1(nkpt1),istwfk2(nkpt2)
 integer,intent(in) :: kg2(3,mpw2*mkmem2),nband1(nkpt1*nsppol1)
 integer,intent(in) :: nband2(nkpt2*nsppol2),ngfft1(18),ngfft2(18)
 integer,intent(in) :: npwarr1(nkpt1),npwarr2(nkpt2),symrel(3,3,nsym)
 real(dp),intent(in) :: gmet1(3,3),gmet2(3,3),kptns1(3,nkpt1),kptns2(3,nkpt2)
 real(dp),intent(in) :: rprimd(3,3),tnons(3,nsym)
 real(dp),intent(inout) :: cg(2,mcg) !vz_i pw_orthon vecnm
 real(dp),intent(inout) :: eigen(mband2*(2*mband2)**formeig*nkpt2*nsppol2)!vz_i newocc
 real(dp),intent(inout) :: occ(mband2*nkpt2*nsppol2) !vz_i

!Local variables-------------------------------
!scalars
 integer,parameter :: init_random=-5,nkpt_max=50,tobox=1,tosph=-1,wr=2
 integer :: aux_stor,band_index,iband,icg,icg_aux,idum
 integer :: ii,ikg2,ikpt1,ikpt10,ikpt2,ikptsp_prev,inplace,iproc
 integer :: isppol1,isppol2,istwf10_k,localrdwf
 integer :: mband1,mband_rd,mband_rw,mcg_aux,me1,me2,mgfft1,mgfft2
 integer :: my_nspinor1,my_nspinor2
 integer :: nb_band,nbd1,nbd1_rd,nbd2,nkpt_eff,nproc2,npw1,npw2,nsp
 integer :: test_cycle,tim_rwwf
 logical :: out_of_core2
 character(len=500) :: message
!arrays
 integer,allocatable :: kg1(:,:),kg2_k(:,:),kg_dum(:,:)
 real(dp) :: kpoint(3),tsec(2)
 real(dp),allocatable :: cg_aux(:,:),eig_k(:),occ_k(:)

! *************************************************************************

 call timab(780,1,tsec)
 call timab(781,1,tsec)

 icg=0

!Init MPI data
 me1=mpi_enreg1%me_kpt
 me2=mpi_enreg2%me_kpt
 nproc2 = mpi_enreg2%nproc_cell
 out_of_core2=(my_nkpt2/=0.and.mkmem2==0)


 if((nsppol1==2.and.nspinor2==2).or.(nspinor1==2.and. nsppol2==2))then
!  This is not yet possible. See later for a message about where to make the needed modifs.
!  EDIT MT 20110707: these modifs are no more needed as they are now done in inwffil
   write(message, '(a,a,a,a,a,a,a,a,a,i2,a,i2,a,a,i2,a,i2,a,a,a,a)' )ch10,&
&   ' newkpt : ERROR -',ch10,&
&   '  The wavefunction translator is (still) unable to interchange',ch10,&
&   '  spin-polarized wfs and spinor wfs. However,',ch10,&
&   '  the input  variables are nsppol1=',nsppol1,', and nspinor1=',nspinor1,ch10,&
&   '  the output variables are nsppol2=',nsppol2,', and nspinor2=',nspinor2,ch10,&
&   '  Action : use a non-spin-polarized wf to start a spinor wf,',ch10,&
&   '           and a non-spinor wf to start a spin-polarized wf.'
   MSG_ERROR(message)
 end if

 my_nspinor1=max(1,nspinor1/mpi_enreg1%nproc_spinor)
 my_nspinor2=max(1,nspinor2/mpi_enreg2%nproc_spinor)
 mband1=maxval(nband1(1:nkpt1*nsppol1))

 if(mkmem1==0 .and. out_of_core2)then
   mband_rd=min(mband1,(mband2/nspinor2)*nspinor1)
   if(mcg<mpw1*my_nspinor1*mband_rd)then
     write(message,'(2(a,i0))')' The dimension mcg= ',mcg,', should be larger than mband_rd= ',mband_rd
     MSG_BUG(message)
   end if
   if(mcg<mband2*mpw2*my_nspinor2)then
     write(message,'(a,i0,a,a,a,i0,a,i0,a,i2)' )&
&     '  The dimension mcg=',mcg,', should be larger than',ch10,&
&     '  the product of mband2=',mband2,', mpw2=',mpw2,', and nspinor2=',my_nspinor2
     MSG_BUG(message)
   end if
 end if

 idum=init_random
 ikpt10 = 0
 istwf10_k=0
 band_index=0
 icg=0

 nkpt_eff=nkpt2
 if( (prtvol==0.or.prtvol==1) .and. nkpt_eff>nkpt_max ) nkpt_eff=nkpt_max

 mgfft1=maxval(ngfft1(1:3))
 mgfft2=maxval(ngfft2(1:3))
 ABI_ALLOCATE(kg1,(3,mpw1))
 ABI_ALLOCATE(kg2_k,(3,mpw2))
 ABI_ALLOCATE(kg_dum,(3,0))

 if (debug>0) then
   if (me1==0) then
     write(std_out,'(a)' ) ' newkpt:  kptns1'
     call prmat (kptns1, 3, nkpt1, 3)
   end if
   if (me2==0) then
     write(std_out,'(a)' ) ' newkpt:  kptns2'
     call prmat (kptns2, 3, nkpt2, 3)
   end if
 end if

 ikptsp_prev=0

 call timab(781,2,tsec)

!Do outer loop over spins
 do isppol2=1,nsppol2

   if (nsppol2==2 .and. me2==0) then
     write(std_out,'(a,i5)' ) ' newkpt: spin channel isppol2 = ',isppol2
   end if

   if (restart==1 .and. out_of_core2) rewind (unkg2)
   ikg2=0

!  Do loop over new k point set
   do ikpt2=1,nkpt2

     call timab(782,1,tsec)

     nbd2=nband2(ikpt2+(isppol2-1)*nkpt2)
     npw2=npwarr2(ikpt2)

     if(restart==1)then

!      Announce the treatment of k point ikpt
       if(ikpt2<=nkpt_eff)then
!        This message might be overwritten in parallel
         write(message, '(a,i6,a,i8,a,i4)' )'P newkpt: treating ',nbd2,' bands with npw=',npw2,' for ikpt=',ikpt2
!        This message might be overwritten in parallel
         if(mpi_enreg2%paralbd==1)then
           do iproc=0,nproc2-1
             nb_band=0
             do iband=1,nbd2
               if(mpi_enreg2%proc_distrb(ikpt2,iband,isppol2) == iproc)nb_band=nb_band+1
             end do
             if(nb_band/=0)then
               write(message, '(a,i6,a,i8,a,i4,a,i4)' ) &
&               'P newkpt: treating ',nb_band,' bands with npw=',npw2,' for ikpt=',ikpt2,' by node ',iproc
             end if
           end do
         end if
         if(mpi_enreg2%paralbd==0) then
           write(message, '(a,i6,a,i8,a,i4,a,i4)' )&
&           'P newkpt: treating ',nbd2,' bands with npw=',npw2,&
&           ' for ikpt=',ikpt2,' by node ',mpi_enreg2%proc_distrb(ikpt2,1,isppol2)
         end if
         if(prtvol>0)then
           call wrtout(iout,message,'COLL')
         end if
       end if

!      Cut the writing if the limit is reached
       if(ikpt2==nkpt_eff+1)then
         if(prtvol>0)then
           call wrtout(iout,' newkpt: prtvol=0 or 1, do not print more k-points.','COLL')
         end if
       end if

!      End of restart==1
     end if

     test_cycle=0
     if(proc_distrb_cycle(mpi_enreg2%proc_distrb,ikpt2,1,nbd2,isppol2,me2)) test_cycle=1
     if(test_cycle==1)then
       if(formeig==0)then
         eigen(1+band_index : nbd2+band_index) = zero
!        occ(1+band_index : nbd2+band_index) = zero
         band_index=band_index+nbd2
       else
         eigen(1+band_index : 2*nbd2**2+band_index) = 0.0_dp
         band_index=band_index+2*nbd2**2
       end if
!      In the case this k point does not belong to me, cycle
       if (my_nkpt2==0) cycle
       if ((mkmem1==0) .and. (ireadwf==1) .and. (mpi_enreg2%paralbd==1))then
         call WffReadSkipK(formeig,headform1,ikpt2,isppol2,mpi_enreg2,wffinp)
         ikptsp_prev=ikptsp_prev+1
       end if
       cycle
     end if

     if(restart==1)then

       if(mkmem2/=0)then
         kg2_k(:,1:npw2)=kg2(:,1+ikg2:npw2+ikg2)
       else if(mkmem2==0)then
!        Read the first line of a block and performs some checks on the unkg file.
         MSG_ERROR("mkmem2 == 0 and rdnpw are not supported anymore.")
         nsp=nspinor2
         !call rdnpw(ikpt2,isppol2,nbd2,npw2,nsp,0,unkg2)
!        Read k+g data
         read (unkg2) kg2_k(1:3,1:npw2)
       end if

     end if

!    Get ikpt1, the closest k from original set, from indkk
     ikpt1=indkk(ikpt2,1)
     if(sppoldbl==2 .and. isppol2==2)ikpt1=indkk(ikpt2+nkpt2,1)

     npw1=npwarr1(ikpt1)
     kpoint(:)=kptns1(:,ikpt1)

!    Determine the spin polarization of the input data
     isppol1=isppol2
     if(nsppol2==2 .and. nsppol1==1)isppol1=1

     if(restart==2)then
       if(ikpt2<=nkpt_eff)then
         write(message,'(a,i4,i8,a,i4,i8)')'- newkpt: read input wf with ikpt,npw=',ikpt1,npw1,', make ikpt,npw=',ikpt2,npw2
         call wrtout(std_out,message,'PERS')
         if(iout/=6 .and. me2==0 .and. prtvol>0)then
           call wrtout(iout,message,'PERS')
         end if
       else if(ikpt2==nkpt_eff+1)then
         write(message,'(a)')'- newkpt : prtvol=0 or 1, do not print more k-points.'
         call wrtout(std_out,message,'PERS')
         if(iout/=6 .and. me2==0 .and. prtvol>0)then
           call wrtout(iout,message,'PERS')
         end if
       end if
     end if

!    Set up the number of bands to be read
     nbd1=nband1(ikpt1+(isppol1-1)*nkpt1)
     nbd1_rd=min(nbd1,(nbd2/nspinor2)*nspinor1)

!    Check that number of bands is not being increased if fill==0 --if so
!    print warning and reset new wf file nband2 to only allowed number
     if ( nbd2/nspinor2 > nbd1/nspinor1 .and. fill==0) then
       if(ikpt2<=nkpt_eff)then
         write(message, '(a,i8,a,i8,a,i8)' )' newkpt: nband2=',nbd2,' < nband1=',nbd1,' => reset nband2 to ',nbd1
         call wrtout(std_out,message,'PERS')
       end if
       nbd2=nbd1
     end if

!    Prepare the reading of the wavefunctions: the correct record is selected
!    WARNING : works only for GS - for RF the number of record differs
     if(restart==2 .and. mkmem1==0)then
       MSG_ERROR("mkmem1 == 0 has been removed.")

       if(debug>0)then
         write(message, '(a,a,a,a,i5,a,i5,a,a,i5,a,i5)' ) ch10,&
         ' newkpt : about to call randac',ch10,&
         '  for ikpt1=',ikpt1,', ikpt2=',ikpt2,ch10,&
         '  and isppol1=',isppol1,', isppol2=',isppol2
         call wrtout(std_out,message,'PERS')
       end if

       !call randac(debug,headform1,ikptsp_prev,ikpt1,isppol1,nband1,nkpt1,nsppol1,wffinp)
     end if

!    Read the data for nbd2 bands at this k point
!    Must decide whether an auxiliary storage is needed
!    When mkmem1==0 and mkmem2==0 , the cg array should be large enough ...
!    When mkmem1==0 and mkmem2/=0 , each k-point block in cg might not be large enough
!    however, will read at most (nbd2/nspinor2)*nspinor1 bands from disk
!    When mkmem1/=0 , it is supposed that each input k-point block is smaller
!    than the corresponding output k-point block, so that the input data
!    have been placed already in cg, at the k-point location where they are needed
     aux_stor=0
     if(mkmem2/=0 .and. mkmem1==0)then
       mcg_aux=npw1*my_nspinor1*nbd1
       if(nbd1_rd<nbd1)mcg_aux=npw1*my_nspinor1*nbd1_rd
       if( mcg_aux > npw2*my_nspinor2*nbd2 )then
         aux_stor=1 ; icg_aux=0
         ABI_ALLOCATE(cg_aux,(2,mcg_aux))
       end if
     end if

     mband_rw=max(nbd1_rd,nbd2)
     ABI_ALLOCATE(eig_k,(mband_rw*(2*mband_rw)**formeig))
     if(formeig==0) then
       ABI_ALLOCATE(occ_k,(mband_rw))
     else
       ABI_ALLOCATE(occ_k,(0))
     end if

     if(mkmem1/=0 .and. ireadwf==1)then
!      Checks that nbd1 and nbd1_rd are equal if eig and occ are input
       if(nbd1/=nbd1_rd)then
         write(message,'(a,a,a,i6,a,i6)')&
&         '  When mkmem1/=0, one must have nbd1=nbd1_rd, while',ch10,&
&         '  nbd1=',nbd1,', and nbd1_rd=',nbd1_rd
         MSG_BUG(message)
       end if
!      Need to put eigenvalues in eig_k, same for occ
!      Note use of band_index, since it is assumed that eigen and occ
!      already have spin-k point location structure than output.
       if(formeig==0)then
         eig_k(1:nbd1_rd)=eigen(1+band_index : nbd1_rd+band_index)
!        occ_k(1:nbd1_rd)=occ(1+band_index : nbd1_rd+band_index)
       else if(formeig==1)then
!        The matrix of eigenvalues has size nbd1 ,  that must be equal
!        to nbd1_rd in the case mkmem1/=0)
         eig_k(1:2*nbd1_rd**2)=eigen(1+band_index : 2*nbd1_rd**2+band_index)
       end if
     end if

     call timab(782,2,tsec)

!    Must read the wavefunctions if they are not yet in place
     if(mkmem1==0 .and. ireadwf==1)then

       if (debug>0 .and. restart==2) then
         write(message,'(a,i5,a,a,i5,a,i5,a)' ) &
&         ' newkpt : about to call rwwf with ikpt1=',ikpt1,ch10,&
&         ' and nband(ikpt1)=',nband1(ikpt1),' nbd2=',nbd2,'.'
         call wrtout(std_out,message,'PERS')
       end if

       if(mpi_enreg1%paralbd==0)tim_rwwf=21
       if(mpi_enreg1%paralbd==1)tim_rwwf=22

       if(aux_stor==0)then
         call rwwf(cg,eig_k,formeig,headform1,icg,ikpt1,isppol1,kg_dum,mband_rw,mcg,mpi_enreg1,&
&         nbd1_rd,nbd1,npw1,my_nspinor1,occ_k,1,0,tim_rwwf,wffinp)
       else
         icg_aux=0
         call rwwf(cg_aux,eig_k,formeig,headform1,icg_aux,ikpt1,isppol1,kg_dum,mband_rw,mcg_aux,&
&         mpi_enreg1,nbd1_rd,nbd1,npw1,my_nspinor1,occ_k,1,0,tim_rwwf,wffinp)
       end if
     end if

     call timab(783,1,tsec)

     if(formeig==1 .and. nbd2/=nbd1_rd .and. ireadwf==1)then
!      Change the storage of eig_k
       if(nbd1_rd<nbd2)then
         do iband=nbd1_rd,1,-1
!          The factor of two is for complex eigenvalues
           do ii=2*nbd2,2*nbd1_rd+1,-1
             eig_k(ii+(iband-1)*2*nbd2)=huge(0.0_dp)/10.0_dp
           end do
           do ii=2*nbd1_rd,1,-1
             eig_k(ii+(iband-1)*2*nbd2)=eig_k(ii+(iband-1)*2*nbd1_rd)
           end do
         end do
       else if(nbd1_rd>nbd2)then
         do iband=1,nbd2
!          The factor of two is for complex eigenvalues
           do ii=1,2*nbd2
             eig_k(ii+(iband-1)*2*nbd2)=eig_k(ii+(iband-1)*2*nbd1_rd)
           end do
         end do
       end if
     end if

!    If change nsppol, must adapt the occupation numbers
!    if(nsppol1/=nsppol2)then
!    occ_k(1:nbd2)=occ_k(1:nbd2)*nsppol1/dbl(nsppol2)
!    then

!    In case nsppol1=2 and nspinor2=2, one should read
!    the other spin-component, and form a spinor wf here, before calling
!    wfconv. One should treat eig_k and occ_k as well.
!    A similar operation is to be performed when nspino1=2 and nsppol2=2
!    EDIT - MT 20110707: the building of the spinor wf is now done in wfffil
!    no need to make it here....

!    DEBUG
!    write(std_out,*)' newkpt: before wfconv'
!    write(std_out,*)' newkpt: mkmem2=',mkmem2
!    stop
!    ENDDEBUG

     call timab(783,2,tsec)
     call timab(784,1,tsec)

!    Note the use of mband2, while mband is used inside
!    write(std_out,*) 'in newkpt,before wfconv,npw1,npw2',npw1,npw2
     inplace=1
     if(aux_stor==0)then
       call wfconv(ceksp2,cg,cg,debug,ecut1,ecut2,ecut2_eff,&
&       eig_k,eig_k,exchn2n3d,formeig,gmet1,gmet2,icg,icg,&
&       ikpt1,ikpt10,ikpt2,indkk,inplace,isppol2,istwfk1,istwfk2,&
&       kg1,kg2_k,kptns1,kptns2,mband_rw,mband_rw,mcg,mcg,&
&       mpi_enreg1,mpi_enreg2,mpw1,mpw2,nbd1_rd,nbd2,&
&       ngfft1,ngfft2,nkpt1,nkpt2,npw1,npw2,nspinor1,nspinor2,nsym,&
&       occ_k,occ_k,optorth,randalg,restart,rprimd,sppoldbl,symrel,tnons)
     else
       call wfconv(ceksp2,cg_aux,cg_aux,debug,ecut1,ecut2,ecut2_eff,&
&       eig_k,eig_k,exchn2n3d,formeig,gmet1,gmet2,icg_aux,icg_aux,&
&       ikpt1,ikpt10,ikpt2,indkk,inplace,isppol2,istwfk1,istwfk2,&
&       kg1,kg2_k,kptns1,kptns2,mband_rw,mband_rw,mcg,mcg,&
&       mpi_enreg1,mpi_enreg2,mpw1,mpw2,nbd1_rd,nbd2,&
&       ngfft1,ngfft2,nkpt1,nkpt2,npw1,npw2,nspinor1,nspinor2,nsym,&
&       occ_k,occ_k,optorth,randalg,restart,rprimd,sppoldbl,symrel,tnons)
     end if

     call timab(784,2,tsec)

!    Finally write new wf to disk file or save in permanent file
     if(mkmem2==0)then

!      Note that in this case, we are sure aux_stor==0
       if(mpi_enreg2%paralbd==0)tim_rwwf=21
       if(mpi_enreg2%paralbd==1)tim_rwwf=22
       call rwwf(cg,eig_k,formeig,0,0,ikpt2,isppol2,kg2_k,nbd2,mcg,mpi_enreg2,&
&       nbd2,nbd2,npw2,my_nspinor2,occ_k,wr,1,tim_rwwf,wffout)

     end if

     call timab(785,1,tsec)

     if(mkmem2/=0)then
       if(aux_stor==1)then
         cg(:,1+icg:npw2*nbd2*my_nspinor2+icg)=cg_aux(:,1:npw2*nbd2*my_nspinor2)
         ABI_DEALLOCATE(cg_aux)
       end if

       icg=icg+npw2*nbd2*my_nspinor2
       ikg2=ikg2+npw2
     end if

     eigen(1+band_index:nbd2*(2*nbd2)**formeig+band_index) = eig_k(1:nbd2*(2*nbd2)**formeig)
!    occ(1+band_index:nbd2+band_index)=occ_k(1:nbd2)

     if(formeig==0)then
       band_index=band_index+nbd2
     else if(formeig==1)then
       band_index=band_index+2*nbd2**2
     end if

     ABI_DEALLOCATE(eig_k)
     ABI_DEALLOCATE(occ_k)

     call timab(785,2,tsec)

   end do ! ikpt2
 end do ! isppol2

 call timab(786,1,tsec)

 if(xmpi_paral==1)then
!  Transmit eigenvalues (not yet occupation numbers)
!  newkpt.F90 is not yet suited for RF format
!  This routine works in both localrdwf=0 or 1 cases.
!  However, in the present routine, localrdwf is to be considered
!  as 1 always, since the transfer has been made in wfsinp .
   localrdwf=1
   call pareigocc(eigen,formeig,localrdwf,mpi_enreg2,mband2,nband2,nkpt2,nsppol2,occ,1)
 end if

 ABI_DEALLOCATE(kg1)
 ABI_DEALLOCATE(kg2_k)
 ABI_DEALLOCATE(kg_dum)

 call timab(786,2,tsec)
 call timab(780,2,tsec)

end subroutine newkpt
!!***
