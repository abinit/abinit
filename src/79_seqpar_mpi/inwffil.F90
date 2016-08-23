!{\src2tex{textfont=tt}}
!!****f* ABINIT/inwffil
!! NAME
!! inwffil
!!
!! FUNCTION
!! Do initialization of wavefunction files.
!! Also call other relevant routines for this initialisation
!! (initialization of wavefunctions from scratch or from file, translations of wavefunctions, ...)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR, AR, MB, MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  ask_accurate= if 1, the wavefunctions and eigenvalues must be
!!    accurate, that is, they must come from a k point that is
!!    symmetric of the needed k point, with a very small tolerance,
!!    the disk file contained sufficient bands to initialize all of them,
!!    the spinor and spin-polarisation characteristics must be identical
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  ecut=effective kinetic energy planewave cutoff (hartree), beyond
!!    which the coefficients of plane waves are zero
!!  ecut_eff=effective kinetic energy planewave cutoff (hartree), needed
!!    to generate the sphere of plane wave
!!  exchn2n3d=if 1, n2 and n3 are exchanged
!!  formeig=explained above
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  ireadwf=option parameter described above for wf initialization
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs to be initialized here.
!!  kg(3,mpw*my_nkpt)=dimensionless coords of G vecs in basis sphere at k point
!!  kptns(3,nkpt)=reduced coords of k points
!!  localrdwf=(for parallel case) if 1, the wffnm  file is local to each machine
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mkmem=number of k-points in core memory
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum number of planewaves as dimensioned in calling routine
!!  nband(nkpt*nsppol)=number of bands at each k point
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt=number of k points
!!  npwarr(nkpt)=array holding npw for each k point.
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of symmetry elements in space group
!!  occ(mband*nkpt*nsppol)=occupations (from disk or left at their initial value)
!!  optorth= 1 if the WFS have to be orthogonalized; 0 otherwise
!!  prtvol=control print volume and debugging
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrel(3,3,nsym)=symmetry operations in real space in terms of primitive translations
!!  tnons(3,nsym)=nonsymmorphic translations for symmetry operations
!!  unkg=unit number for storage of basis sphere data: stores indirect
!!   indexing array and integer coordinates for all planewaves in basis
!!   sphere for each k point being considered
!!  unwff1,unwfnow= unit numbers for files wffnm and wft1nm.
!!  wffnm=name (character data) of file for input wavefunctions.
!!
!! OUTPUT
!!  wff1  = structure information for files wffnm .
!!  wffnow= structure information for wf file wft1nm
!!  if ground state format (formeig=0):
!!    eigen(mband*nkpt*nsppol)=eigenvalues (input or init to large number), (Ha)
!!  if respfn format (formeig=1):
!!    eigen(2*mband*mband*nkpt*nsppol)=matrix of eigenvalues
!!                                     (input or init to large number), (Ha)
!! Conditional output (returned if mkmem/=0):
!!  cg(2,mcg)=complex wf array
!!    be careful : an array of size cg(2,npw*nspinor), as used
!!    in the response function code, is not enough !
!!  wvl <type(wvl_data)>=all wavelets data.
!!
!! NOTES
!! Detailed description :
!!  Initialize unit wff1%unwff for input of wf data if ireadwf=1
!!  Opens file on unit wffnow%unwff
!!   if the storage on disk is needed (mkmem==0)
!!  Initializes wf data on wffnow%unwff, by calling the appropriate routine.
!!
!! formeig option (format of the eigenvalues and occupations) :
!!   0 => ground-state format (initialisation of
!!        eigenvectors with random numbers, vector of eigenvalues,
!!        occupations are present)
!!   1 => respfn format (initialisation of
!!        eigenvectors with 0 s, hermitian matrix of eigenvalues)
!!
!! ireadwf options:
!!   0 => initialize with random numbers or 0 s
!!   1 => read from disk file wff1, initializing higher bands
!!        with random numbers or 0 s if not provided in disk file
!!
!! The wavefunctions after this initialisation are stored in unit wffnow%unwff
!!
!! WARNINGS
!! The symmetry operations are used to translate the data from one
!! k point to another, symmetric, k point.
!! They can be completely different from the symmetry operations
!! contained on the disk file. No check is performed between the two sets.
!!
!! Occupations will not be modified nor output, in the present status of this routine.
!!
!! If ground state format (formeig=0) occ(mband*nkpt*nsppol) was output.
!! NOT OUTPUT NOW !
!!
!! PARENTS
!!      dfpt_looppert,dfptnl_loop,gstate,nonlinear,respfn
!!
!! CHILDREN
!!      copy_mpi_enreg,destroy_mpi_enreg,hdr_check,hdr_free,hdr_io,hdr_ncread
!!      kpgio,listkk,matr3inv,newkpt,pawrhoij_copy,timab,wffopen,wfsinp,wrtout
!!      wvl_wfsinp_disk,wvl_wfsinp_scratch
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine inwffil(ask_accurate,cg,dtset,ecut,ecut_eff,eigen,exchn2n3d,&
&           formeig,gmet,hdr,ireadwf,istwfk,kg,kptns,localrdwf,mband,&
&           mcg,mkmem,mpi_enreg,mpw,nband,ngfft,nkpt,npwarr,&
&           nsppol,nsym,occ,optorth,rprimd,symafm,symrel,tnons,unkg,wff1,&
&           wffnow,unwff1,wffnm,wvl)

 use defs_basis
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_wffile
 use m_errors
 use m_xmpi
 use m_nctk
 use m_hdr

 use m_io_tools, only : file_exists
 use m_pawrhoij, only : pawrhoij_type, pawrhoij_copy, pawrhoij_io
 use m_mpinfo,   only : destroy_mpi_enreg, copy_mpi_enreg

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'inwffil'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_56_recipspace
 use interfaces_67_common
 use interfaces_79_seqpar_mpi, except_this_one => inwffil
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: ask_accurate,exchn2n3d,formeig,ireadwf,localrdwf,mband,mcg,mkmem,mpw
 integer,intent(in) :: nkpt,nsppol,nsym,optorth,unkg,unwff1
 real(dp),intent(in) :: ecut,ecut_eff
 character(len=*),intent(in) :: wffnm
 type(MPI_type),intent(inout),target :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(wffile_type),intent(inout) :: wff1
 type(wffile_type),intent(inout) :: wffnow !vz_i
 type(wvl_data),intent(inout) :: wvl
 integer,intent(in) :: istwfk(nkpt),kg(3,mpw*mkmem),ngfft(18)
 integer,intent(in) :: npwarr(nkpt),symafm(nsym),symrel(3,3,nsym)
 integer,intent(in),target :: nband(nkpt*nsppol)
 real(dp),intent(inout),target :: cg(2,mcg),eigen((2*mband)**formeig*mband*nkpt*nsppol) !vz_i
 real(dp),intent(in) :: gmet(3,3),kptns(3,nkpt),rprimd(3,3),tnons(3,nsym)
 real(dp),intent(inout) :: occ(mband*nkpt*nsppol)

!Local variables-------------------------------
 integer,parameter :: master=0
 integer :: iomode,accurate,ceksp,debug,doorth,fform,fform_dum,fill
 integer :: headform0,iband,ibg,ibg0,icg,icg0,icgsft,ieigsft,ierr,ii
 integer :: ikassoc,ikpt,ikpt0,ikptsp,ikptsp0,imax,increase_nkassoc,isppol,isppol0
 integer :: mband0,mband0_rd,mband_eff,mcg_disk,me,me0,mkmem0,mpw0
 integer :: my_nkpt,my_nspinor,my_nspinor0,nband_k,nband0_k
 integer :: nkassoc,nkpt0,npw,npw0,nspinor0,nspinor_eff,nsppol0,nsppol_eff,nsppol2nspinor
 integer :: rdwr,randalg,restart,restartpaw,spaceComm,spaceComm_io,sppoldbl,sppoldbl_eff,squeeze
 logical :: out_of_core
 real(dp) :: dksqmax,ecut0
 character(len=500) :: message
 type(hdr_type) :: hdr0
 integer :: ngfft0(18)
 integer,allocatable :: indkk0(:,:),indx(:),istwfk0(:),kg0(:,:)
 integer,allocatable :: nband0_rd(:),npwarr0(:),npwi(:),npwtot0(:)
 integer,allocatable,target :: indkk(:,:),nband0(:)
 integer, pointer :: indkk_eff(:,:),nband_eff(:)
 logical,allocatable :: my_kpt(:)
 real(dp) :: gmet0(3,3),gprim0(3,3),rprim0(3,3),tsec(2)
 real(dp),allocatable :: cg_disk(:,:),kptns0(:,:)
 real(dp),pointer :: cg_eff(:,:),eigen_eff(:)
 type(MPI_type),pointer :: mpi_enreg0

! *************************************************************************

 DBG_ENTER("COLL")

!Keep track of total time spent in inwffil
 call timab(710,1,tsec)
 call timab(711,1,tsec)

!Check the validity of formeig
 if (formeig/=0.and.formeig/=1) then
   write(message,'(a,i0,a)')' formeig = ',formeig,', but the only allowed values are 0 or 1.'
   MSG_BUG(message)
 end if

!Init mpi_comm
 spaceComm=mpi_enreg%comm_cell
 spaceComm_io=xmpi_comm_self
 if (mpi_enreg%paral_kgb==1) spaceComm_io= mpi_enreg%comm_bandspinorfft
 if (mpi_enreg%paral_hf ==1) spaceComm_io= mpi_enreg%comm_hf
 me=xmpi_comm_rank(spaceComm)

!Determine number of k points processed by current node
 my_nkpt=nkpt;if (size(mpi_enreg%my_kpttab)>0) my_nkpt=maxval(mpi_enreg%my_kpttab)
 out_of_core=(mkmem==0.and.my_nkpt/=0)

 ngfft0(:)=ngfft(:)
 headform0=0 !Default value for headform0 (will be needed later, to read wf blocks)

 ! Chebyshev is more sensitive to the quality of input random numbers, so use a new algorithm
 if(dtset%wfoptalg == 1) then
   randalg = 1
 else
   ! Otherwise, use compatibility mode
   randalg = 0
 end if

!If the input data are on disk, determine the kind of restart
 wff1%fname = wffnm

 !  Checking the existence of data file
 if (ireadwf==1 .and. .not.file_exists(wff1%fname)) then
   ! Trick needed to run Abinit test suite in netcdf mode.
   if (file_exists(nctk_ncify(wff1%fname))) then
     write(std_out,"(3a)")"- File: ",trim(wff1%fname)," does not exist but found netcdf file with similar name."
     wff1%fname = nctk_ncify(wff1%fname)
   end if
   if (localrdwf/=0 .and. .not. file_exists(wff1%fname)) then
     MSG_ERROR('Missing data file: '//TRIM(wff1%fname))
   end if
 end if

 if (ireadwf==1)then

   iomode=dtset%iomode
   if (localrdwf==0) then
     ! This is in case the wff file must be read by only the master proc
     if (iomode /= IO_MODE_ETSF) iomode=IO_MODE_FORTRAN_MASTER 
     !iomode=IO_MODE_FORTRAN_MASTER 
   end if

   call WffOpen(iomode,spaceComm,wff1%fname,ierr,wff1,master,me,unwff1,spaceComm_io)

!  Initialize hdr0 (sent to all procs), thanks to reading of wff1
   rdwr=1
   if ( ANY(wff1%iomode == (/IO_MODE_FORTRAN_MASTER, IO_MODE_FORTRAN, IO_MODE_MPI/) )) then
     call hdr_io(fform_dum,hdr0,rdwr,wff1)
#ifdef HAVE_TRIO_NETCDF
   else if (wff1%iomode == IO_MODE_ETSF) then
     call hdr_ncread(hdr0, wff1%unwff, fform_dum)
#endif
   end if

   ! Handle IO Error.
   if (fform_dum == 0) then
     write(message,"(4a)")&
&     "hdr_io returned fform == 0 while trying to read the wavefunctions from file: ",trim(wff1%fname),ch10,&
&     "This usually means that the file does not exist or that you don't have enough privileges to read it"
     MSG_ERROR(message)
   end if

   call wrtout(std_out,'inwffil: examining the header of disk file '//TRIM(wff1%fname),'COLL')

!  Check hdr0 versus hdr (and from now on ignore header consistency and write new info to header for each file)
   if (dtset%usewvl == 0) then
!    wait for plane waves.
     fform=2
   else
!    wait for wavelets.
     fform = 200
   end if
   call hdr_check(fform,fform_dum,hdr,hdr0,'PERS',restart,restartpaw)

   nkpt0=hdr0%nkpt
   nsppol0=hdr0%nsppol
   headform0=hdr0%headform

   write(message,'(2a)')'-inwffil : will read wavefunctions from disk file ',trim(wff1%fname)
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

 else
   restart=1; restartpaw=0

!  Fill some data concerning an hypothetical file to be read
!  This is to allow the safe use of same routines than with ireadwf==1.
   nkpt0=nkpt ; nsppol0=nsppol
 end if ! end ireadwf

 sppoldbl=1
 if(minval(symafm(:))==-1)then
   if(nsppol0==1 .and. nsppol==2)sppoldbl=2
 end if

 ABI_ALLOCATE(indkk,(nkpt*sppoldbl,6))
 ABI_ALLOCATE(istwfk0,(nkpt0))
 ABI_ALLOCATE(kptns0,(3,nkpt0))
 ABI_ALLOCATE(nband0,(nkpt0*nsppol0))
 ABI_ALLOCATE(npwarr0,(nkpt0))

 if(restart==2)then ! restart with translations

   ecut0=hdr0%ecut_eff
   istwfk0(1:nkpt0)=hdr0%istwfk(1:nkpt0)
   kptns0(1:3,1:nkpt0)=hdr0%kptns(1:3,1:nkpt0)
   nband0(1:nkpt0*nsppol0)=hdr0%nband(1:nkpt0*nsppol0)
   ngfft0(1:3)=hdr0%ngfft(1:3)
   npwarr0(1:nkpt0)=hdr0%npwarr(1:nkpt0)
   nspinor0=hdr0%nspinor
   rprim0(:,:)=hdr0%rprimd(:,:)
   mpw0=maxval(npwarr0(:))

!  Compute reciprocal space metric gmet for unit cell of disk wf
   call matr3inv(rprim0,gprim0)
   do ii=1,3
     gmet0(:,ii)=gprim0(1,:)*gprim0(1,ii)+gprim0(2,:)*gprim0(2,ii)+gprim0(3,:)*gprim0(3,ii)
   end do

   if ((mpi_enreg%paral_kgb==1).or.(mpi_enreg%paral_hf==1)) then
     ABI_DATATYPE_ALLOCATE(mpi_enreg0,)
     call copy_mpi_enreg(mpi_enreg,mpi_enreg0)
     ABI_ALLOCATE(kg0,(3,mpw0*nkpt0))
     ABI_ALLOCATE(npwtot0,(nkpt0))
     message="tmpfil"
     call kpgio(ecut0,dtset%exchn2n3d,gmet0,istwfk0,kg0, &
&     kptns0,nkpt0,nband0,nkpt0,'PERS',mpi_enreg0,&
&     mpw0,npwarr0,npwtot0,nsppol0)

     ABI_DEALLOCATE(kg0)
     ABI_DEALLOCATE(npwtot0)
   else
     mpi_enreg0 => mpi_enreg
   end if

!  At this stage, the header of the file wff1i%unwff is read, and
!  the pointer is ready to read the first wavefunction block.

!  Compute k points from input file closest to the output file
   call listkk(dksqmax,gmet0,indkk,kptns0,kptns,nkpt0,nkpt,nsym,sppoldbl,symafm,symrel,1)

 else if (restart==1) then ! direct restart

!  Fill variables that must be the same, as determined by hdr_check.f
!  This is to allow the safe use of the same routines than with restart==2.
   nspinor0=dtset%nspinor
   ecut0=ecut_eff
   gmet0(:,:)=gmet(:,:)
   istwfk0(:)=istwfk(:)
   kptns0(:,:)=kptns(:,:)
   npwarr0(:)=npwarr(:)
   mpw0=mpw

   do isppol=1,sppoldbl
     do ikpt=1,nkpt
       indkk(ikpt+(isppol-1)*nkpt,1)=ikpt
       indkk(ikpt+(isppol-1)*nkpt,2:6)=0
     end do
   end do
   dksqmax=0.0_dp

!  The treatment of nband0 asks for some care
   if(ireadwf==0)then
     nband0(:)=0
   else
     nband0(1:nkpt0*nsppol0)=hdr0%nband(1:nkpt0*nsppol0)
   end if

   mpi_enreg0 => mpi_enreg

 else
   mpi_enreg0 => mpi_enreg
 end if

 if(mpi_enreg0%paral_pert == 1.and.mpi_enreg0%me_pert/=-1) then
   me0 = mpi_enreg0%me_pert
 else
   me0 = mpi_enreg0%me_cell
 end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Before hdr_free:
!If restartpaw==1, store hdr0%pawrhoij in hdr%pawrhoij; else if restartpaw==0,
!hdr%pawrhoij(:)has been initialized in hdr_init.
 if(restartpaw==1) then
   call pawrhoij_copy(hdr0%pawrhoij,hdr%pawrhoij,keep_itypat=.true.)
 end if

 call timab(711,2,tsec)
 call timab(712,1,tsec)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!At this stage, all the relevant information from the header of the disk file,
!has been exploited, and stored in variables, on all processors.
!It is also contained in hdr0
!(on all processors, except if restart=1 and localrdwf=0,
!in which case it is only on the master)
!These information might be changed later, while processing the
!wavefunction data, and converting it. The variable hdr0 might be kept
!for further checking, or reference, or debugging, but at present,
!it is simpler to close it. The other header, hdr, will be used for the new file, if any.

 if(ask_accurate==1)then

!  Check whether the accuracy requirements might be fulfilled
   if(ireadwf==0)then
     write(message,'(9a)')&
&     'The file ',trim(wff1%fname),' cannot be used to start the ',ch10,&
&     'present calculation. It was asked that the wavefunctions be accurate,',ch10,&
&     'but they were not even read.',ch10,&
&     'Action: use a wf file, with ireadwf/=0.'
     MSG_ERROR(message)
   end if
   if(dksqmax>tol12)then
     write(message, '(9a,es16.6,4a)' )&
&     'The file ',trim(wff1%fname),' cannot be used to start the ',ch10,&
&     'present calculation. It was asked that the wavefunctions be accurate, but',ch10,&
&     'at least one of the k points could not be generated from a symmetrical one.',ch10,&
&     'dksqmax=',dksqmax,ch10,&
&     'Action: check your wf file and k point input variables',ch10,&
&     '        (e.g. kptopt or shiftk might be wrong in the present dataset or the preparatory one.'
     MSG_ERROR(message)
   end if
   if(dtset%nspinor/=nspinor0)then
     write(message,'(a,a, a,a,a,a,a, a,a,2i5,a,a)')&
&     'The file ',trim(wff1%fname),' cannot be used to start the ',ch10,&
&     'present calculation. It was asked that the wavefunctions be accurate, but',ch10,&
&     'nspinor differs in the file from the actual nspinor.',ch10,&
&     'nspinor,nspinor0=',dtset%nspinor,nspinor0,ch10,&
&     'Action: check your wf file, and nspinor input variables.'
     MSG_ERROR(message)
   end if
   if((nsppol>nsppol0 .and. sppoldbl==1) .or. nsppol<nsppol0 ) then
     write(message,'(a,a, a,a,a,a,a, a,a,3i5,a,a)')&
&     'The file ',trim(wff1%fname),' cannot be used to start the ',ch10,&
&     'present calculation. It was asked that the wavefunctions be accurate, but',ch10,&
&     'the nsppol variables do not match in the file and in the actual calculation',ch10,&
&     'nsppol,nsppol,sppoldbl=',dtset%nspinor,nspinor0,sppoldbl,ch10,&
&     'Action: check your wf file, and nsppol input variables.'
     MSG_ERROR(message)
   end if

!  Now, check the number of bands
   accurate=1
   do isppol=1,nsppol
     do ikpt=1,nkpt
       ikpt0=indkk(ikpt+(isppol-1)*(sppoldbl-1)*nkpt,1)
       ikptsp =ikpt +(isppol-1)*nkpt
       ikptsp0=ikpt0+(isppol-1)*(2-sppoldbl)*nkpt0
       if(nband0(ikptsp0)<nband(ikptsp))accurate=0
     end do
   end do
   if(accurate==0)then
     write(message,'(a,a, a,a,a,a,a, a,a)')&
&     'The file ',trim(wff1%fname),' cannot be used to start the ',ch10,&
&     'present calculation. It was asked that the wavefunctions be accurate,',ch10,&
&     'but the number of bands differ in the file and in the actual calculation.',ch10,&
&     'Action: use a wf file with the correct characteristics.'
     MSG_ERROR(message)
   end if

 end if

!Flag: do we need to translate WF to (from) spinors ?
 nsppol2nspinor=0
 if (nsppol0==2.and.dtset%nspinor==2) nsppol2nspinor=+1
 if (nspinor0==2.and.nsppol==2) nsppol2nspinor=-1

!Take into account parallism over spinors
 my_nspinor =max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 my_nspinor0=max(1,nspinor0/mpi_enreg0%nproc_spinor)

!Not all bands might be read, if not needed to fill the wavefunctions
 mband0=maxval(nband0(1:nkpt0*nsppol0))
 mband0_rd=min(mband0,(mband/dtset%nspinor)*nspinor0)

!****************************************************************************
!If needed, transfer the input wf from disk to core memory
!(in the parallel case, it allows to change localrdwf=0 in localrdwf=1)

 mkmem0=0

 if(xmpi_paral == 1 .or. mpi_enreg%paral_kgb == 1 .or. mpi_enreg%paral_hf == 1) then
   if(localrdwf==0 .and. out_of_core)then
     MSG_BUG('localrdwf==0 and mkmem==0 (out-of-core solution) are not allowed together (yet)')
   end if
 end if

 call timab(712,2,tsec)

!Here, treat reading wavefunctions with mkmem/=0, first step
 if(ireadwf==1 .and. (.not.out_of_core))then

   call timab(713,1,tsec)

!  if(restart==1 .and. ireadwf==1 .and. mkmem/=0)then

!  Compute table of k point associations. Make a trial choice for nkassoc.
   nkassoc=(nkpt/nkpt0+1)*2
   ABI_ALLOCATE(indkk0,(nkpt0,nkassoc))
!  Infinite loops are allowed in F90
   do
     indkk0(:,:)=0
     increase_nkassoc=0
     do ikpt=1,nkpt*sppoldbl
       ikpt0=indkk(ikpt,1)
       do ikassoc=1,nkassoc
         if(indkk0(ikpt0,ikassoc)==0)then
           indkk0(ikpt0,ikassoc)=ikpt
           exit
         end if
         if(nkassoc==ikassoc)increase_nkassoc=1
       end do
       if(increase_nkassoc==1)then
         ABI_DEALLOCATE(indkk0)
         nkassoc=2*nkassoc
         ABI_ALLOCATE(indkk0,(nkpt0,nkassoc))
         exit
       end if
     end do
     if(increase_nkassoc==0)exit
   end do

!  DEBUG
!  write(std_out,*)' inwffil: indkk0, nkassoc=',nkassoc
!  do ikpt0=1,nkpt0
!  write(std_out,*)' ikpt0,indkk0(ikpt0,1)=',ikpt0,indkk0(ikpt0,1)
!  end do
!  ENDDEBUG

!  DEBUG
!  write(std_out,*)' inwffil : indkk(:,1)=',indkk(:,1)
!  write(std_out,*)' inwffil : sppoldbl=',sppoldbl
!  ENDDEBUG

!  To treat the case (nsppol0=2,nspinor0=1)<->(nsppol=1,nspinor=2),
!  apply the following trick:
!  1- We call wfsinp with fake arguments (nsppol_eff and nspinor_eff)
!  2- We transform collinear polarized WF into spinors
!  or  spinors into collinear polarized WF
   if (nsppol2nspinor/=0.and.out_of_core.and.dtset%usewvl==0) then
     write(message, '(7a)')&
&     'When mkmem=0 (out-of-core), the wavefunction translator is unable',ch10,&
&     'to interchange spin-polarized wfs and spinor wfs.',ch10,&
&     'Action: use a non-spin-polarized wf to start a spinor wf,',ch10,&
&     '        and a non-spinor wf to start a spin-polarized wf.'
     MSG_ERROR(message)
   end if

!  === Fake arguments definition for wfsinp
   if (nsppol2nspinor==0.or.dtset%usewvl/=0) then
     indkk_eff => indkk
     nband_eff => nband
     eigen_eff => eigen
     cg_eff => cg
     nspinor_eff=dtset%nspinor;nsppol_eff=nsppol;sppoldbl_eff=sppoldbl
     mband_eff=maxval(nband_eff(1:nkpt*nsppol_eff))
   else if (nsppol2nspinor==1.and.(.not.out_of_core)) then
     nsppol_eff=2;nspinor_eff=1;sppoldbl_eff=1
     ABI_ALLOCATE(indkk_eff,(nkpt*sppoldbl_eff,6))
     ABI_ALLOCATE(nband_eff,(nkpt*nsppol_eff))
     indkk_eff(1:nkpt,1:6)   =indkk(1:nkpt,1:6)
     nband_eff(1:nkpt)       =nband(1:nkpt)/2
     nband_eff(1+nkpt:2*nkpt)=nband(1:nkpt)/2
     mband_eff=maxval(nband_eff(1:nkpt*nsppol_eff))
     eigen_eff => eigen
     cg_eff => cg
   else if (nsppol2nspinor==-1.and.(.not.out_of_core)) then
!    WARNING: MT 07072011 -> this is memory consuming
!    A copy a spinorial WF (and eigenvalues) is temporary kept in memory;
!    But the case (nspinor=2 => nsppol=2) might be rare
!    and only useful for testing purposes.
!    => print a warning for the user
!    NOTE: in that case (nsppol=2), parallelization over spinors is not activated

     write(message,'(5a)')&
&     'In the case of spinor WF read from disk and converted into',ch10,&
&     'spin-polarized non-spinor WF, the WF translator is memory',ch10,&
&     'consuming (a copy a spinor WF is temporary stored in memory).'
     MSG_WARNING(message)

     nsppol_eff=1;nspinor_eff=2;sppoldbl_eff=1
     ABI_ALLOCATE(indkk_eff,(nkpt*sppoldbl_eff,6))
     ABI_ALLOCATE(nband_eff,(nkpt*nsppol_eff))
     indkk_eff(1:nkpt,1:6)=indkk(1:nkpt,1:6)
     nband_eff(1:nkpt)    =2*nband(1:nkpt)
     mband_eff=maxval(nband_eff(1:nkpt*nsppol_eff))
     ABI_ALLOCATE(eigen_eff,((2*mband_eff)**formeig*mband_eff*nkpt*nsppol_eff))
     ABI_ALLOCATE(cg_eff,(2,mpw0*nspinor_eff*mband_eff*mkmem*nsppol_eff))
   end if

!  === nband0 argument definition for wfsinp
   squeeze=0
   ABI_ALLOCATE(cg_disk,(0,0))
   if(.not.out_of_core)then
     ABI_ALLOCATE(nband0_rd,(nkpt0*nsppol0))
     nband0_rd(:)=0
     do isppol=1,nsppol_eff
       do ikpt=1,nkpt
         ikpt0=indkk(ikpt+(isppol-1)*(sppoldbl-1)*nkpt,1)
         isppol0=min(isppol,nsppol0)
         ikptsp =ikpt +(isppol -1)*nkpt
         ikptsp0=ikpt0+(isppol0-1)*(2-sppoldbl)*nkpt0
         nband0_k=min(nband0(ikptsp0),(nband_eff(ikptsp)/nspinor_eff)*nspinor0)
         nband0_rd(ikptsp0)=max(nband0_rd(ikptsp0),nband0_k)
         npw0=npwarr0(ikpt0)
         npw =npwarr (ikpt)
         if(npw0*nspinor0*nband0_k > npw*nspinor_eff*nband_eff(ikptsp))squeeze=1
       end do
     end do
     if(squeeze==1)then
       mcg_disk=mpw0*my_nspinor0*mband0_rd
       ABI_DEALLOCATE(cg_disk)
       ABI_ALLOCATE(cg_disk,(2,mcg_disk))
     else
       if(xmpi_paral == 1 .or. mpi_enreg0%paral_kgb == 1 .or. mpi_enreg0%paral_hf == 1)then
         if(localrdwf==0)then
           mcg_disk=mpw0*my_nspinor0*mband0_rd
           ABI_DEALLOCATE(cg_disk)
           ABI_ALLOCATE(cg_disk,(2,mcg_disk))
         end if
       end if
     end if
   end if

   call timab(713,2,tsec)
   call timab(714,1,tsec)

!  === call to wfsinp
   if (dtset%usewvl == 0) then
     call wfsinp(cg_eff,cg_disk,ecut,ecut0,ecut_eff,eigen,&
&     exchn2n3d,formeig,gmet,gmet0,headform0,&
&     indkk_eff,indkk0,istwfk,istwfk0,kptns,kptns0,localrdwf,&
&     mband_eff,mcg,mcg_disk,mpi_enreg,mpi_enreg0,mpw,mpw0,&
&     nband_eff,nband0_rd,ngfft,nkassoc,nkpt,nkpt0,npwarr,npwarr0,nspinor_eff,nspinor0,&
&     nsppol_eff,nsppol0,nsym,occ,optorth,dtset%prtvol,randalg,restart,rprimd,sppoldbl_eff,squeeze,&
&     symrel,tnons,wff1)
     if (nsppol2nspinor/=0)  then
       ABI_DEALLOCATE(indkk_eff)
       ABI_DEALLOCATE(nband_eff)
     end if
   else
!    Read wavefunctions from file.
     call wvl_wfsinp_disk(dtset, hdr0, hdr, mpi_enreg, occ, 1, &
&     hdr%rprimd, wff1, wvl%wfs, wvl%descr, hdr%xred)
   end if

   call timab(714,2,tsec)
   call timab(715,1,tsec)

!  Now, update xyz0 variables, for use in newkpt
   nband0(:)=nband0_rd(:)

!  If squeeze, the conversion was done in wfsinp, so no conversion left.
   if(squeeze==1)then
     ecut0=ecut_eff
     gmet0(:,:)=gmet(:,:)
     ABI_DEALLOCATE(kptns0)
     ABI_DEALLOCATE(istwfk0)
     ABI_DEALLOCATE(nband0)
     ABI_DEALLOCATE(npwarr0)
     ABI_ALLOCATE(kptns0,(3,nkpt))
     ABI_ALLOCATE(istwfk0,(nkpt))
     ABI_ALLOCATE(nband0,(nkpt*nsppol))
     ABI_ALLOCATE(npwarr0,(nkpt))
     kptns0(:,:)=kptns(:,:)
     istwfk0(:)=istwfk(:)
     npwarr0(:)=npwarr(:)
     nband0(:)=0
     do isppol=1,nsppol
       do ikpt=1,nkpt
         ikpt0=indkk(ikpt+(isppol-1)*(sppoldbl-1)*nkpt,1)
         isppol0=min(isppol,nsppol0)
         ikptsp =ikpt +(isppol -1)*nkpt
         ikptsp0=ikpt0+(isppol0-1)*(sppoldbl-1)*nkpt0
         nband0(ikptsp)=(nband0_rd(ikptsp0)/nspinor0)*dtset%nspinor
       end do
     end do
     do ikpt=1,nkpt
       indkk(ikpt,1)=ikpt
       indkk(ikpt,2:6)=0
     end do
!    This transfer must come after the nband0 transfer
     nspinor0=dtset%nspinor
     nkpt0=nkpt
     nsppol0=nsppol
   end if ! end squeeze == 1

!  The input wavefunctions have been transferred from disk to core memory
   mkmem0=mkmem

   ABI_DEALLOCATE(indkk0)
   ABI_DEALLOCATE(nband0_rd)
   ABI_DEALLOCATE(cg_disk)

   call timab(715,2,tsec)

 else !ireadwf == 0
   if (dtset%usewvl == 1) then

     call timab(714,1,tsec)
!    Compute wavefunctions from input guess.
     call wvl_wfsinp_scratch(dtset, mpi_enreg, occ, hdr%rprimd, wvl, hdr%xred)
     call timab(714,2,tsec)
   end if
 end if

 call timab(716,1,tsec)

!=== Eventual conversion of WF into (from) spinors
 if (dtset%usewvl==0) then

!  ***** No conversion (standard case) ****
   if (nsppol2nspinor==0) then
     nspinor_eff=nspinor0;nsppol_eff=nsppol0;sppoldbl_eff=sppoldbl
     indkk_eff => indkk
     nband_eff => nband0

!    ***** Conversion from collinear to spinorial WF ****
   else if (nsppol2nspinor==1.and.(.not.out_of_core)) then
!    Translate the WF and eigenvalues from nsppol=2 to nspinor=2
!    This is tricky (because we do not want to create a temporary array for cg)
     nsppol_eff=1;nspinor_eff=2;sppoldbl_eff=1
     ABI_ALLOCATE(indkk_eff,(nkpt*sppoldbl_eff,6))
     ABI_ALLOCATE(nband_eff,(nkpt0*nsppol_eff))
     indkk_eff(1:nkpt,1:6)=indkk(1:nkpt,1:6)
     nband_eff(1:nkpt0)=2*nband0(1:nkpt0)
!    Compute some shifts from isspol0=1 to isppol0=2
     imax=0;icgsft=0;ieigsft=0
     ABI_ALLOCATE(my_kpt,(nkpt0))
     do ikpt0=1,nkpt0
       nband0_k=nband0(ikpt0);nband_k=nband(ikpt0)
       my_kpt(ikpt0)=(.not.(proc_distrb_cycle(mpi_enreg0%proc_distrb,ikpt0,1,nband_k,1,me0)))
       ieigsft=ieigsft+(2*nband0_k)**formeig*nband0_k
       if(my_kpt(ikpt0)) then
         imax=imax+nband0_k;icgsft=icgsft+nband0_k*npwarr0(ikpt0)
       end if
     end do
!    --- First version: no parallelization over spinors
     if (mpi_enreg0%paral_spinor==0) then
!      Compute some useful indexes
       ABI_ALLOCATE(indx,(2*imax))
       ABI_ALLOCATE(npwi,(imax))
       ii=0;icg=0
       do ikpt0=1,nkpt0
         if(my_kpt(ikpt0)) then
           nband0_k=nband0(ikpt0);npw0=npwarr0(ikpt0)
           do iband=1,nband0_k
             ii=ii+1;npwi(ii)=npw0
             indx(2*ii-1)=icg+mpw0;indx(2*ii)=icg+2*mpw0
             icg=icg+4*mpw0
           end do
         end if
       end do
!      Expand WF in cg (try to use the whole array)
       ii=nsppol0*imax;icg0=nsppol0*icgsft
       do isppol=nsppol0,1,-1
         do ikpt0=nkpt0,1,-1
           if(my_kpt(ikpt0)) then
             nband0_k=nband0(ikpt0);npw0=npwarr0(ikpt0)
             do iband=nband0_k,1,-1
               icg0=icg0-npw0 
               if (indx(ii)<icg0) then 
                 MSG_BUG("Unable to read WF!")
               end if
               cg(:,indx(ii)+1:indx(ii)+npw0)=cg(:,icg0+1:icg0+npw0)
               ii=ii-1
             end do
           end if
         end do
       end do
!      Convert polarized WF into spinors
       ii=1
       do ikpt0=1,nkpt0
         if(my_kpt(ikpt0)) then
           nband0_k=nband0(ikpt0);npw0=npwarr0(ikpt0)
           do iband=1,nband0_k
             npw0=npwi(ii)
             cg(:,indx(2*ii-1)-mpw0+1:indx(2*ii-1)-mpw0+npw0)=cg(:,indx(ii)+1:indx(ii)+npw0)
             cg(:,indx(2*ii  )+mpw0+1:indx(2*ii  )+mpw0+npw0)=cg(:,indx(ii+imax)+1:indx(ii+imax)+npw0)
             ii=ii+1
           end do
         end if
       end do
!      Compress new cg array (from mpw to npw) and cancel zero-components
       icg0=0;icg=0
       do ikpt0=1,nkpt0
         if(my_kpt(ikpt0)) then
           nband0_k=nband0(ikpt0);npw0=npwarr0(ikpt0)
           do iband=1,nband0_k
             cg(:,icg0       +1:icg0+  npw0)=cg(:,icg+1:icg+npw0)
             cg(:,icg0+  npw0+1:icg0+2*npw0)=zero
             cg(:,icg0+2*npw0+1:icg0+3*npw0)=zero
             cg(:,icg0+3*npw0+1:icg0+4*npw0)=cg(:,icg+3*mpw0+1:icg+3*mpw0+npw0)
             icg0=icg0+4*npw0;icg=icg+4*mpw0
           end do
         end if
       end do
!      --- Second version: parallelization over spinors
     else
!      Compute some useful indexes
       ABI_ALLOCATE(indx,(imax))
       ABI_ALLOCATE(npwi,(imax))
       ii=0;icg=0
       do ikpt0=1,nkpt0
         if(my_kpt(ikpt0)) then
           nband0_k=nband0(ikpt0);npw0=npwarr0(ikpt0)
           do iband=1,nband0_k
             ii=ii+1;npwi(ii)=npw0
             indx(ii)=icg+mpi_enreg0%me_spinor*mpw0
             icg=icg+2*mpw0
           end do
         end if
       end do
!      Expand WF in cg
       ii=(mpi_enreg0%me_spinor+1)*imax;icg0=(mpi_enreg0%me_spinor+1)*icgsft
       do ikpt0=nkpt0,1,-1
         if(my_kpt(ikpt0)) then
           nband0_k=nband0(ikpt0);npw0=npwarr0(ikpt0)
           do iband=nband0_k,1,-1
             icg0=icg0-npw0 
             if (indx(ii)<icg0) then 
               MSG_BUG("Unable to read WF!")
             end if
             cg(:,indx(ii)+1:indx(ii)+npw0)=cg(:,icg0+1:icg0+npw0)
             ii=ii-1
           end do
         end if
       end do
!      Compress new cg array (from mpw to npw) and cancel zero-components
       icg0=0;icg=0
       do ikpt0=1,nkpt0
         if(my_kpt(ikpt0)) then
           nband0_k=nband0(ikpt0);npw0=npwarr0(ikpt0)
           do iband=1,nband0_k
             if (mpi_enreg0%me_spinor==0) then
               cg(:,icg0     +1:icg0+  npw0)=cg(:,icg+1:icg+npw0)
               cg(:,icg0+npw0+1:icg0+2*npw0)=zero
             else
               cg(:,icg0     +1:icg0+  npw0)=zero
               cg(:,icg0+npw0+1:icg0+2*npw0)=cg(:,icg+mpw0+1:icg+mpw0+npw0)
             end if
             icg0=icg0+2*npw0;icg=icg+2*mpw0
           end do
         end if
       end do
     end if
!    Translate eigenvalues
     ibg0=2*ieigsft;ibg=2*ieigsft
     do ikpt0=nkpt0,1,-1
       nband0_k=nband0(ikpt0)
       ibg0=ibg0-  nband0_k*(2*nband0_k)**formeig
       ibg =ibg -2*nband0_k*(2*nband0_k)**formeig
       if(my_kpt(ikpt0)) then
         do iband=nband0_k*(2*nband0_k)**formeig,1,-1
           eigen(2*iband-1+ibg)=eigen(iband+ibg0-ieigsft)
           eigen(2*iband  +ibg)=eigen(iband+ibg0)
         end do
       end if
     end do
     ABI_DEALLOCATE(indx)
     ABI_DEALLOCATE(npwi)
     ABI_DEALLOCATE(my_kpt)

!    ***** Conversion from spinorial to collinear WF ****
   else if (nsppol2nspinor==-1.and.(.not.out_of_core)) then
!    In that case parallelization over spinors is never activated
     nsppol_eff=2;nspinor_eff=1;sppoldbl_eff=1
     ABI_ALLOCATE(indkk_eff,(nkpt*sppoldbl_eff,6))
     ABI_ALLOCATE(nband_eff,(nkpt0*nsppol_eff))
     indkk_eff(1:nkpt,1:6)=indkk(1:nkpt,1:6)
     nband_eff(1:nkpt0)        =nband0(1:nkpt0)/2
     nband_eff(1+nkpt0:2*nkpt0)=nband0(1:nkpt0)/2
!    Compute shifts from isspol0=1 to isppol0=2
     icgsft=0;ieigsft=0
     do ikpt0=1,nkpt0
       nband0_k=nband0(ikpt0);nband_k=nband(ikpt0)
       ieigsft=ieigsft+(nband0_k/2)*(nband0_k)**formeig
       if(.not.(proc_distrb_cycle(mpi_enreg0%proc_distrb,ikpt0,1,nband_k,1,me))) &
&       icgsft=icgsft+(nband0_k/2)*npwarr0(ikpt0)
     end do
!    Translate the WF and eigenvalues from nspinor=2 to nsppol=2
     icg0=0;icg=0;ibg=0
     do ikpt0=1,nkpt0
       nband0_k=nband0(ikpt0);nband_k=nband(ikpt0);npw0=npwarr0(ikpt0)
       if(.not.(proc_distrb_cycle(mpi_enreg0%proc_distrb,ikpt0,1,nband_k,1,me))) then
         do iband=1,nband0_k/2
           do ii=1,npw0
             cg(:,ii+icg)       =cg_eff(:,ii+icg0)
             cg(:,ii+icg+icgsft)=cg_eff(:,ii+icg0+3*npw0)
           end do
           icg0=icg0+4*npw0;icg=icg+npw0
         end do
         do iband=(nband0_k/2)*(nband0_k)**formeig,1,-1
           eigen(iband+ibg)        =eigen_eff(2*iband-1+2*ibg)
           eigen(iband+ibg+ieigsft)=eigen_eff(2*iband  +2*ibg)
!          occ(iband+ibg)        =occ_eff(2*iband-1+2*ibg)
!          occ(iband+ibg+ieigsft)=occ_eff(2*iband  +2*ibg)
         end do
       end if
       ibg=ibg+(nband0_k/2)*(nband0_k)**formeig
     end do
     ABI_DEALLOCATE(cg_eff)
     ABI_DEALLOCATE(eigen_eff)

   else
     MSG_BUG('unable to interchange nsppol and nspinor when mkmem=0')
   end if
 end if

!Clean hdr0
 !if (ireadwf==1)then
 !  if( restart==2 .or. localrdwf==1 .or. master==me)then
 !    call hdr_free(hdr0)
 !  end if
 !end if
 call hdr_free(hdr0)

 call timab(716,2,tsec)
 call timab(717,1,tsec)


!****************************************************************************
!Now, treat translation of wavefunctions if wavefunctions are planewaves

 ceksp=0; debug=0; doorth=1; fill=1
 if (dtset%usewvl == 0) then
   
   call newkpt(ceksp,cg,debug,ecut0,ecut,ecut_eff,eigen,exchn2n3d,&
&   fill,formeig,gmet0,gmet,headform0,indkk_eff,&
&   ab_out,ireadwf,istwfk0,istwfk,kg,kptns0,kptns,&
&   mband,mcg,mkmem0,mkmem,mpi_enreg0,mpi_enreg,&
&   mpw0,mpw,my_nkpt,nband_eff,nband,ngfft0,ngfft,nkpt0,nkpt,npwarr0,npwarr,&
&   nspinor_eff,dtset%nspinor,nsppol_eff,nsppol,nsym,occ,optorth,&
&   dtset%prtvol,randalg,restart,rprimd,sppoldbl_eff,symrel,tnons,unkg,wff1,wffnow)

   if (nsppol2nspinor/=0)  then
     ABI_DEALLOCATE(nband_eff)
   end if

 end if ! dtset%usewvl == 0

!****************************************************************************

 ABI_DEALLOCATE(indkk)
 ABI_DEALLOCATE(istwfk0)
 ABI_DEALLOCATE(kptns0)
 ABI_DEALLOCATE(nband0)
 ABI_DEALLOCATE(npwarr0)
 if (restart==2 .and.(mpi_enreg0%paral_kgb==1 .or. mpi_enreg0%paral_hf == 1)) then
   call destroy_mpi_enreg(mpi_enreg0)
   ABI_DATATYPE_DEALLOCATE(mpi_enreg0)
 else
   nullify(mpi_enreg0)
 end if

 call timab(717,2,tsec)
 call timab(710,2,tsec)

 DBG_EXIT("COLL")

end subroutine inwffil
!!***
