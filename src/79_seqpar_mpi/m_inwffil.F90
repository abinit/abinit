!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_inwffil
!! NAME
!!  m_inwffil
!!
!! FUNCTION
!!  Do initialization of wavefunction files.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA, XG, GMR, AR, MB, MVer, ZL, MB, TD)
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

module m_inwffil

 use defs_basis
 use defs_wvltypes
 use m_abicore
 use m_wffile
 use m_wfk
 use m_errors
 use m_xmpi
 use m_nctk
 use m_hdr
 use m_dtset
#if defined HAVE_MPI2
 use mpi
#endif

 use defs_abitypes, only : MPI_type
 use m_time,     only : timab
 use m_io_tools, only : file_exists, get_unit
 use m_geometry, only : getspinrot
 use m_pptools,  only : prmat
 use m_symtk,    only : matr3inv, mati3inv
 use m_cgtools,  only : cg_envlop, pw_orthon
 use m_fftcore,  only : kpgsph, sphere, sphereboundary
 use m_pawrhoij, only : pawrhoij_type, pawrhoij_copy, pawrhoij_io
 use m_mpinfo,   only : destroy_mpi_enreg, copy_mpi_enreg, proc_distrb_cycle
 use m_kg,       only : kpgio, ph1d3d, getph
 use m_kpts,     only : listkk
 use m_occ,      only : pareigocc
 use m_rwwf,     only : rwwf, WffReadSkipK
 use m_wvl_wfsinp, only : wvl_wfsinp_disk, wvl_wfsinp_scratch

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

 private
!!***

 public ::  inwffil
!!***

contains
!!***

!!****f* m_inwffil/inwffil
!! NAME
!! inwffil
!!
!! FUNCTION
!! Do initialization of wavefunction files.
!! Also call other relevant routines for this initialisation
!! (initialization of wavefunctions from scratch or from file, translations of wavefunctions, ...)
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
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  ireadwf=option parameter described above for wf initialization
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs to be initialized here.
!!  kg(3,mpw*my_nkpt)=dimensionless coords of G vecs in basis sphere at k point
!!  kptns(3,nkpt)=reduced coords of k points
!!  localrdwf=(for parallel case) if 1, the wffnm  file is local to each machine
!!  mband=maximum number of bands
!!  mband_mem=maximum number of bands for this cpu
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband_mem*mkmem*nsppol
!!  mkmem=number of k-points in core memory
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum number of planewaves as dimensioned in calling routine
!!  nband(nkpt*nsppol)=number of bands at each k point
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkpt=number of k points
!!  npwarr(nkpt)=array holding npw for each k point.
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of symmetry elements in space group
!!  occ(mband*nkpt*nsppol)=occupations (from disk or left at their initial value)
!!  optorth= 1 if the WFS have to be orthogonalized; 0 otherwise
!!  prtvol=control print volume and debugging
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

subroutine inwffil(ask_accurate,cg,dtset,ecut,ecut_eff,eigen,exchn2n3d,&
&           formeig,hdr,ireadwf,istwfk,kg,kptns,localrdwf,mband,&
&           mcg,mkmem,mpi_enreg,mpw,nband,ngfft,nkpt,npwarr,&
&           nsppol,nsym,occ,optorth,symafm,symrel,tnons,unkg,wff1,&
&           wffnow,unwff1,wffnm,wvl)

!Arguments ------------------------------------
 integer,intent(in) :: ask_accurate,exchn2n3d,formeig,ireadwf,localrdwf,mband,mcg,mkmem,mpw
 integer,intent(in) :: nkpt,nsppol,nsym,optorth,unkg,unwff1
 real(dp),intent(in) :: ecut,ecut_eff
 character(len=*),intent(in) :: wffnm
 type(MPI_type),intent(inout),target :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(wffile_type),intent(inout) :: wff1
 type(wffile_type),intent(inout) :: wffnow
 type(wvl_data),intent(inout) :: wvl
 integer,intent(in) :: istwfk(nkpt),kg(3,mpw*mkmem),ngfft(18)
 integer,intent(in) :: npwarr(nkpt),symafm(nsym),symrel(3,3,nsym)
 integer,intent(in),target :: nband(nkpt*nsppol)
 real(dp),intent(inout),target :: cg(2,mcg),eigen((2*mband)**formeig*mband*nkpt*nsppol)
 real(dp),intent(in) :: kptns(3,nkpt),tnons(3,nsym)
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
 real(dp) :: gmet(3,3),gmet0(3,3),gprim0(3,3),rprim0(3,3),tsec(2)
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

!Chebyshev is more sensitive to the quality of input random numbers, so use a new algorithm
 if(dtset%wfoptalg == 1) then
   randalg = 1
 else
   ! Otherwise, use compatibility mode
   randalg = 0
 end if

!If the input data are on disk, determine the kind of restart
 wff1%fname = wffnm

!Checking the existence of data file
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

!Compute reciprocal space metric gmet
 call matr3inv(hdr%rprimd,gprim0) ! gprim0 is used as temporary storage
 gmet=matmul(transpose(gprim0),gprim0)

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
#ifdef HAVE_NETCDF
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

   call wrtout(std_out,' inwffil: examining the header of disk file: '//trim(wff1%fname),'COLL')

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
   gmet0=matmul(transpose(gprim0),gprim0)

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
   call listkk(dksqmax,gmet0,indkk,kptns0,kptns,nkpt0,nkpt,nsym,sppoldbl,symafm,symrel,1,spaceComm)

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
&     'consuming (a copy of the spinor WF is temporarily stored in memory).'
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
&     nsppol_eff,nsppol0,nsym,occ,optorth,dtset%prtvol,randalg,restart,hdr%rprimd,sppoldbl_eff,squeeze,&
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
 call hdr0%free()

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
&   dtset%prtvol,randalg,restart,hdr%rprimd,sppoldbl_eff,symrel,tnons,unkg,wff1,wffnow)

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

!!****f* m_inwffil/wfsinp
!! NAME
!! wfsinp
!!
!! FUNCTION
!! Do initialization of wavefunction files.
!! Also call other relevant routines for this initialisation.
!! Detailed description :
!!  - Initialize unit wff1 for input of wf data
!!
!! formeig option (format of the eigenvalues and occupations) :
!!   0 => ground-state format (initialisation of
!!        eigenvectors with random numbers, vector of eigenvalues,
!!        occupations are present)
!!   1 => respfn format (initialisation of
!!        eigenvectors with 0 s, hermitian matrix of eigenvalues)
!!
!! INPUTS
!!  ecut0=kinetic energy cutoffs for basis sphere 0 (hartree) (if squeeze=1)
!!  ecut=kinetic energy cutoffs beyond which the coefficients of cg vanish (Ha)
!!   (needed only if squeeze=1)
!!  ecut_eff=effective kinetic energy planewave cutoff (hartree), needed
!!    to generate the sphere of plane wave
!!  exchn2n3d=if 1, n2 and n3 are exchanged
!!  formeig=explained above
!!  gmet(3,3), gmet0(3,3)=reciprocal space metrics (bohr^-2)
!!  headform0=header format (might be needed to read the block of wfs)
!!  indkk(nkpt*sppoldbl,6)=describe k point number of kptns0 that allows to
!!   generate wavefunctions closest to given kpt
!!   indkk(:,1)=k point number of kptns0
!!   indkk(:,2)=symmetry operation to be applied to kpt0, to give kpt0a
!!    (if 0, means no symmetry operation, equivalent to identity )
!!   indkk(:,3:5)=shift in reciprocal space to be given to kpt0a,
!!    to give kpt0b, that is the closest to kpt.
!!   indkk(:,6)=1 if time-reversal was used to generate kpt1a from kpt1, 0 otherwise
!!  indkk0(nkpt0,nkassoc)=list of k points that will be generated by k point number ikpt0
!!  istwfk(nkpt)=input parameter that describes the storage of wfs
!!  istwfk0(nkpt0)=input parameter that describes the storage of wfs in set0
!!  kptns(3,nkpt),kptns0(3,nkpt0)=k point sets (reduced coordinates)
!!  localrdwf=(for parallel case) if 1, the wff1%unwff file is local to each machine
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*nsppol
!!  mpi_enreg=information about MPI parallelization
!!  mpi_enreg0=information about MPI parallelization in set0
!!  mpw=maximum number of planewaves as dimensioned in calling routine
!!  mpw0=maximum number of planewaves on disk file
!!  nban_dp_rd(nkpt0*nsppol0)=number of bands to be read at each k point
!!  nband(nkpt*nsppol)=number of bands at each k point
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkassoc=dimension of indkk0 array
!!  nkpt=number of k points expected
!!  nkpt0=number of k points on disk
!!  npwarr(nkpt)=array holding npw for each k point.
!!  npwarr0(nkpt0)=array holding npw for each k point, disk format.
!!  nspinor=number of spinorial components of the wavefunctions
!!  nspinor0=number of spinorial components of the wavefunctions on disk
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsppol0=1 for unpolarized, 2 for spin-polarized, when on disk
!!  nsym=number of symmetry elements in space group
!!  optorth= 1 if the WFS have to be orthogonalized; 0 otherwise
!!  prtvol=control print volume and debugging
!!  randalg=1 if "good" (but non-portable) random numbers should be used, 0 for compatibility
!!  restart= if 2, want conversion between wavefunctions
!!           if 1, direct restart is allowed (see hdr_check.f)
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  sppoldbl= if 1, no doubling of the number if spins thanks to antiferromagn
!!            if 2, deduce nsppol=2 from nsppol=1, using Shubnikov symmetries
!!  squeeze=1 if cg_disk is to be used, even with mkmem/=0
!!  symrel(3,3,nsym)=symmetry operations in real space in terms
!!   of primitive translations
!!  tnons(3,nsym)=nonsymmorphic translations for symmetry operations
!!  wff1, structure information for input and output files
!!
!! OUTPUT
!!  if ground state format (formeig=0):
!!    eigen(mband*nkpt*nsppol)=eigenvalues (input or init to large number), (Ha)
!!  if respfn format (formeig=1):
!!    eigen(2*mband*mband*nkpt*nsppol)=
!!         matrix of eigenvalues (input or init to large number), (Ha)
!! Conditional output:
!!    cg_disk(2,mpw*nspinor*mband*mkmem*nsppol)=complex wf array
!!      be careful : an array of size cg(2,npw*nspinor), as used
!!      in the response function code, is not enough !
!!
!! SIDE EFFECTS
!!  if ground state format (formeig=0):
!!    occ(mband*nkpt*nsppol)=occupations (from disk or left at their initial value)
!!    NOT OUTPUT NOW !
!!
!! NOTES
!! occ will not be modified nor output, in the present status of this routine.
!!
!! WARNINGS
!! For parallelism : no distinction yet between nban_dp_rd and nband
!!
!! TODO
!! THE DESCRIPTION IS TO BE COMPLETELY REVISED, AS THIS ONE COMES FROM inwffil.f
!!
!! PARENTS
!!      inwffil
!!
!! CHILDREN
!!      initwf,mpi_bcast,mpi_recv,mpi_send,pareigocc,timab,wfconv,wffreadskipk
!!      wrtout
!!
!! SOURCE

subroutine wfsinp(cg,cg_disk,ecut,ecut0,ecut_eff,eigen,exchn2n3d,&
&                  formeig,gmet,gmet0,headform0,indkk,indkk0,istwfk,&
&                  istwfk0,kptns,kptns0,localrdwf,mband,&
&                  mcg,mcg_disk,mpi_enreg,mpi_enreg0,mpw,mpw0,nband,nban_dp_rd,&
&                  ngfft,nkassoc,nkpt,nkpt0,npwarr,npwarr0,nspinor,&
&                  nspinor0,nsppol,nsppol0,nsym,occ,optorth,prtvol,randalg,restart,rprimd,&
&                  sppoldbl,squeeze,symrel,tnons,wff1)

!Arguments ------------------------------------
 integer, intent(in) :: exchn2n3d,formeig,headform0,localrdwf,mband,mcg,mcg_disk
 integer, intent(in) :: mpw,mpw0,nkassoc,nkpt,nkpt0,nspinor,nspinor0,nsppol,nsppol0,nsym
 integer, intent(in) :: optorth,prtvol,randalg,restart,sppoldbl,squeeze
 real(dp), intent(in) :: ecut,ecut0,ecut_eff
 type(MPI_type), intent(inout) :: mpi_enreg,mpi_enreg0
 type(wffile_type), intent(inout) :: wff1
 integer, intent(in) :: indkk(nkpt*sppoldbl,6),indkk0(nkpt0,nkassoc),istwfk(nkpt)
 integer, intent(in) :: istwfk0(nkpt0),nband(nkpt*nsppol),nban_dp_rd(nkpt0*nsppol0)
 integer, intent(in) :: ngfft(18),npwarr(nkpt),npwarr0(nkpt0),symrel(3,3,nsym)
 real(dp), intent(in) :: gmet(3,3),gmet0(3,3),kptns(3,nkpt),kptns0(3,nkpt0),rprimd(3,3)
 real(dp), intent(in) :: tnons(3,nsym)
 real(dp), intent(out) :: eigen((2*mband)**formeig*mband*nkpt*nsppol)
 real(dp), intent(inout) :: cg(2,mcg),cg_disk(2,mcg_disk) !vz_i pw_ortho
 real(dp), intent(inout) :: occ(mband*nkpt*nsppol)

!Local variables-------------------------------
 integer :: band_index,band_index_trial,ceksp,debug,dim_eig_k,iband,icg
 integer :: icg_disk,icg_trial,idum,ierr,ii,ikassoc,ikassoc_trial,ikpt,ikpt0
 integer :: ikpt10,ikpt_trial,ikptsp,ikptsp_old,inplace,isp,isp_max,isppol,isppol0
! integer :: ipw ! commented out below
 integer :: isppol_trial,me,mgfft,my_nspinor,my_nspinor0
 integer :: nban_dp_k,nban_dp_rdk,nband_k,nband_rdk,nband_trial,nbd,nbd_max
 integer :: ncopy,nkpt_eff,nproc_max,npw0_k,npw_k,npw_ktrial
 integer :: read_cg,read_cg_disk,sender,spaceComm
 character(len=500) :: message
 integer,allocatable :: band_index_k(:,:),icg_k(:,:),kg0_k(:,:),kg_k(:,:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: eig0_k(:),eig_k(:),occ0_k(:),occ_k(:)
#if defined HAVE_MPI
 integer :: iproc,my_ikpt
 integer :: tag,test_cycle
 integer :: statux(MPI_STATUS_SIZE)
 integer,allocatable :: ikassoc_me(:),ikpt_me(:),isppol_me(:),nband_k_me(:)
#endif
 integer :: nkpt_max=50

! *************************************************************************

!DEBUG
!write(std_out,*)' wfsinp : enter'
!write(std_out,*)' wfsinp : nband=',nband(:)
!write(std_out,*)' wfsinp : nban_dp_rd=',nban_dp_rd(:)
!write(std_out,*)' wfsinp : localrdwf=',localrdwf
!write(std_out,*)' wfsinp : paralbd,formeig=',mpi_enreg%paralbd,formeig
!write(std_out,*)' wfsinp : indkk0(:,1)=',indkk0(:,1)
!ENDDEBUG

 call timab(720,1,tsec)
 call timab(721,3,tsec)

 nkpt_max=50; if(xmpi_paral==1)nkpt_max=-1
 nbd_max=size(mpi_enreg%proc_distrb,2)
 isp_max=size(mpi_enreg%proc_distrb,3)

!Init mpi_comm
 spaceComm=mpi_enreg%comm_cell
 nproc_max=xmpi_comm_size(spaceComm)
 me=mpi_enreg%me_kpt
 sender = 0

#if defined HAVE_MPI
 if(localrdwf==0)then
   ABI_ALLOCATE(ikpt_me,(nproc_max))
   ABI_ALLOCATE(nband_k_me,(nproc_max))
   ABI_ALLOCATE(ikassoc_me,(nproc_max))
   ABI_ALLOCATE(isppol_me,(nproc_max))
 end if
#endif

!Check the validity of formeig
 if(formeig/=0.and.formeig/=1)then
   write(message, '(a,i0,a)' )' formeig=',formeig,' , but the only allowed values are 0 or 1.'
   MSG_BUG(message)
 end if

 my_nspinor =max(1,nspinor /mpi_enreg%nproc_spinor)
 my_nspinor0=max(1,nspinor0/mpi_enreg%nproc_spinor)

 nkpt_eff=max(nkpt0,nkpt)
 if( (prtvol==0.or.prtvol==1) .and. nkpt_eff>nkpt_max)nkpt_eff=nkpt_max

 ABI_ALLOCATE(icg_k,(nkpt,nsppol))
 ABI_ALLOCATE(band_index_k,(nkpt,nsppol))

!write(std_out,*)' wfsinp : me,isppol,ikpt,icg_k,band_index_k'

!Compute the locations of the blocks in cg, eig and occ
 icg=0
 band_index=0
 ikpt10=0

 do isppol=1,nsppol
   do ikpt=1,nkpt

     nband_k=nband(ikpt+(isppol-1)*nkpt)
     band_index_k(ikpt,isppol)=band_index

#if defined HAVE_MPI
     test_cycle=0;nbd=min(nband_k,nbd_max);isp=min(isppol,isp_max)
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nbd,isp,me))test_cycle=1
     if(test_cycle==1)then
       band_index=band_index+nband_k*(2*nband_k)**formeig
!      In the case this k point does not belong to me, cycle
       cycle
     end if
#endif

     npw_k=npwarr(ikpt)
     icg_k(ikpt,isppol)=icg
     icg=icg+npw_k*my_nspinor*nband_k

     band_index=band_index+nband_k*(2*nband_k)**formeig
!    write(std_out,'(5i8)' )me,isppol,ikpt,icg_k(ikpt,isppol),band_index_k(ikpt,isppol)
   end do ! End k point loop
 end do  ! End spin loop

 band_index=0
 ikptsp_old=0

!DEBUG
!write(std_out,*)' wfsinp: before loop'
!write(std_out,*)' nsppol0,nsppol,nkpt0',nsppol0,nsppol,nkpt0
!write(std_out,*)' mpw,mgfft,mpw,mpw0',mpw,mgfft,mpw,mpw0
!ENDDEBUG

 mgfft=maxval(ngfft(1:3))
 if(squeeze==1)then
   ABI_ALLOCATE(kg_k,(3,mpw))
   ABI_ALLOCATE(kg0_k,(3,mpw0))
 end if

 eigen(:)=0.0_dp
!occ(:)=0.0_dp

 call timab(721,2,tsec)

!Loop over spins
!For the time being, do not allow nsppol=2 to nspinor=2 conversion
!MT 20110707: this can be done by a fake call to the routine: see inwffil
 do isppol0=1,min(nsppol0,nsppol)

!  Loop on k points  : get the cg then eventually write on unwfnow
   do ikpt0=1,nkpt0

     call timab(722,1,tsec)

     nban_dp_rdk=nban_dp_rd(ikpt0+(isppol0-1)*nkpt0)

!    DEBUG
!    write(std_out,*)' wfsinp: ikpt0,isppol0,nkpt0=',ikpt0,isppol0,nkpt0
!    write(std_out,*)' nban_dp_rdk=',nban_dp_rdk
!    ENDDEBUG

     npw0_k=npwarr0(ikpt0)
     if(ikpt0<=nkpt_eff)then
       write(message,'(a,a,2i4)')ch10,' wfsinp: inside loop, init ikpt0,isppol0=',ikpt0,isppol0
       call wrtout(std_out,message,'PERS')
     end if

!    Must know whether this k point is needed, and in which
!    block (ikpt, isppol), the wavefunction is to be placed.
!    Select the one for which the number of bands is the biggest.
     ikpt=0
     isppol=0
     ikassoc=0
     nband_k=0
#if defined HAVE_MPI
     if(localrdwf==0)then
       nband_k_me(:)=0
       ikpt_me(:)=0
       isppol_me(:)=0
       ikassoc_me(:)=0
       nband_k_me(:)=0
     end if
#endif

     do isppol_trial=1,nsppol

       if(nsppol==2 .and. nsppol0==2 .and. isppol0/=isppol_trial)cycle

       do ikassoc_trial=1,nkassoc

         ikpt_trial=indkk0(ikpt0,ikassoc_trial)
         if(sppoldbl==2)then
           if(isppol_trial==1 .and. ikpt_trial>nkpt)cycle
           if(isppol_trial==2 .and. ikpt_trial<=nkpt)cycle
           if(isppol_trial==2)ikpt_trial=ikpt_trial-nkpt
         end if

#if defined HAVE_MPI
         if(localrdwf==1)then
           if(ikpt_trial/=0)then
             nband_trial=nband(ikpt_trial+(isppol_trial-1)*nkpt)
             nbd=min(nband_trial,nbd_max);isp=min(isppol_trial,isp_max)
             if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt_trial,1,nbd,isp,me))ikpt_trial=0
           end if
         end if
#endif

         if(ikpt_trial/=0)then
           nband_trial=nband(ikpt_trial+(isppol_trial-1)*nkpt)
           if(nband_k<nband_trial)then
             nband_k=nband_trial ; ikpt=ikpt_trial ; isppol=isppol_trial
             ikassoc=ikassoc_trial
           end if

#if defined HAVE_MPI
           if(localrdwf==0)then
             do iproc=1,nproc_max
               my_ikpt=1
               nband_trial=nband(ikpt_trial+(isppol_trial-1)*nkpt)
               nbd=min(nband_trial,nbd_max);isp=min(isppol_trial,isp_max)
               if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt_trial,1,nbd,isp,(iproc-1))) my_ikpt=0
               if(my_ikpt/=0)then
                 nband_trial=nband(ikpt_trial+(isppol_trial-1)*nkpt)
                 if(nband_k_me(iproc)<nband_trial)then
                   nband_k_me(iproc)=nband_trial
                   ikpt_me(iproc)=ikpt_trial
                   isppol_me(iproc)=isppol_trial
                   ikassoc_me(iproc)=ikassoc_trial
                 end if
               end if
             end do
           end if
#endif

         end if

       end do ! ikassoc_trial
     end do ! isppol_trial

!    DEBUG
!    write(std_out,*)' wfsinp : me,select isppol,ikpt=',me,isppol,ikpt
#if defined HAVE_MPI
!    write(std_out,*)' wfsinp : me,ikpt_me(:)=',me,ikpt_me(:)
#endif
!    write(std_out,*)' wfsinp : me,isppol_me(:)=',me,isppol_me(:)
!    stop
!    ENDDEBUG

     call timab(722,2,tsec)

!    If the wavefunction block to be read is interesting ...
     if (ikpt/=0)then

       call timab(723,3,tsec)
       sender = me
       nband_k=nband(ikpt+(isppol-1)*nkpt)
       npw_k=npwarr(ikpt)

#if defined HAVE_MPI
       if (localrdwf==1.or.(localrdwf==0.and.me==0)) then

         if(ikpt<=nkpt_eff)then
           write(message,'(a,i6,a,i8,a,i4,a,i4)') &
&           ' wfsinp: treating ',nband_k,' bands with npw=',npw_k,' for ikpt=',ikpt,' by node ',me
           call wrtout(std_out,message,'PERS')
         else if(ikpt==nkpt_eff+1)then
           call wrtout(std_out,' wfsinp: prtvol=0 or 1, do not print more k-points.','PERS')
         end if

       end if
#endif

       nband_rdk=nban_dp_rdk
       if(squeeze==1)nband_rdk=(nban_dp_rdk/nspinor0)*nspinor
       if(formeig==0)then
         ABI_ALLOCATE(eig_k,(nband_rdk))
         ABI_ALLOCATE(occ_k,(nband_rdk))
         ABI_ALLOCATE(eig0_k,(nban_dp_rdk))
         ABI_ALLOCATE(occ0_k,(nban_dp_rdk))
         dim_eig_k=nband_rdk
       else if(formeig==1)then
         ABI_ALLOCATE(eig_k,(2*nband_rdk*nband_rdk))
         ABI_ALLOCATE(eig0_k,(2*nban_dp_rdk*nban_dp_rdk))
         dim_eig_k=2*nband_rdk*nband_rdk
         ABI_ALLOCATE(occ0_k,(0))
         ABI_ALLOCATE(occ_k,(0))
       else
         ABI_ALLOCATE(occ0_k,(0))
         ABI_ALLOCATE(eig0_k,(0))
         ABI_ALLOCATE(occ_k,(0))
         ABI_ALLOCATE(eig_k,(0))
       end if
       eig_k(:)=0.0_dp
       eig0_k(:)=0.0_dp

!      Generate or read the cg for this k point
!      Either read into cg, or read into cg_disk
       read_cg=1 ; read_cg_disk=0
       if(squeeze==1)then
         read_cg=0 ; read_cg_disk=1
       end if
#if defined HAVE_MPI
       if(localrdwf==0)then
         read_cg=0
         read_cg_disk=0
!        XG20040106 The following condition is correct
         if(me==0)read_cg_disk=1
       end if
#endif

       icg=0
       if(read_cg==1)icg=icg_k(ikpt,isppol)

!      DEBUG
!      write(std_out,*)' wfsinp: before initwf',wff1%offwff
!      write(std_out,*)' wfsinp: me,read_cg,read_cg_disk=',me,read_cg,read_cg_disk
!      write(std_out,*)' wfsinp: nban_dp_rdk=',nban_dp_rdk
!      ENDDEBUG
       call timab(723,2,tsec)
       call timab(724,3,tsec)

       if(read_cg_disk==1)then
         call initwf (cg_disk,eig0_k,formeig,headform0,icg,ikpt0,ikptsp_old,&
&         isppol0,mcg_disk,mpi_enreg0, &
&         nban_dp_rdk,nkpt0,npw0_k,my_nspinor0,occ0_k,wff1)
       end if

       if(read_cg==1)then
         call initwf (cg,eig0_k,formeig,headform0,icg,ikpt0,ikptsp_old,&
&         isppol0,mcg,mpi_enreg0,&
&         nban_dp_rdk,nkpt0,npw0_k,my_nspinor0,occ0_k,wff1)
       end if

       call timab(724,2,tsec)
       call timab(725,3,tsec)

       nban_dp_k=min(nban_dp_rdk,(nband_k/nspinor)*nspinor0)
!      This band_index is defined BEFORE the eventual redefinition
!      of ikpt and isppol, needed  when localrdwf==0 in parallel
       band_index=band_index_k(ikpt,isppol)
!      DEBUG
!      write(std_out,*)' wfsinp: me,cg_disk(:,1)=',me,cg_disk(:,1)
!      ENDDEBUG

!      DEBUG
!      if(me==0 .and. ikpt0==1)then
!      write(std_out,*)' wfsinp : cg array, before trial, ikpt0=',ikpt0
!      do ipw=1,15
!      write(std_out,'(i4,2es20.10)' )ipw,cg(:,ipw)
!      end do
!      end if
!      ENDDEBUG


#if defined HAVE_MPI
       if(localrdwf==0)then
!        Warning: In that case , not yet // on nspinors
!        Transmit to each of the other processors, when needed
         if(nproc_max>=2)then
           do iproc=2,nproc_max
!            Only me=0 and me=iproc-1 are concerned by this
             if(me==0 .or. me==iproc-1)then

               ikpt=ikpt_me(iproc)
               isppol=isppol_me(iproc)

               if(ikpt/=0)then
!                In this case, processor iproc-1 needs the data
!                Generate a common tag
                 tag=256*(ikpt-1)+iproc+1
                 if(isppol==2)tag=-tag
                 nband_k=nband(ikpt+(isppol-1)*nkpt)
                 npw_k=npwarr(ikpt)
!                SEND
                 if(me==0)then
                   write(std_out,*)'SENDWFSINP ',me
                   call MPI_SEND(cg_disk,2*npw_k*my_nspinor*nband_k,&
&                   MPI_DOUBLE_PRECISION,iproc-1,tag,spaceComm,ierr)
                 end if
!                RECEIVE
                 if(me==iproc-1)then
                   call MPI_RECV(cg_disk,2*npw_k*my_nspinor*nband_k,&
&                   MPI_DOUBLE_PRECISION,0,tag,spaceComm,statux,ierr)
                   icg=icg_k(ikpt,isppol)
                   if(squeeze==0)then
                     cg(:,icg+1:icg+npw_k*my_nspinor*nband_k)=&
&                     cg_disk(:,1:npw_k*my_nspinor*nband_k)
                   end if
                   ikassoc=ikassoc_me(iproc)
                 end if
               end if

             end if
           end do ! iproc
         end if

!        DEBUG
!        write(std_out,*)' wfsinp: me, iproc loop finished',me
!        ENDDEBUG

!        Take care of me=0 needing the data
         if (me==0) then
           ikpt=ikpt_me(me+1)
           isppol=isppol_me(me+1)
           if(ikpt/=0 )then
             nband_k=nband(ikpt+(isppol-1)*nkpt)
             npw_k=npwarr(ikpt)
!            I am the master node, and I might need my own data
             icg=icg_k(ikpt,isppol)
             if(squeeze==0)then
!              Copy from cg_disk to cg
               cg(:,1+icg:npw_k*my_nspinor*nband_k+icg)= &
&               cg_disk(:,1:npw_k*my_nspinor*nband_k)
             end if
             ikassoc=ikassoc_me(me+1)
           end if
         end if
!        For the eigenvalues and occ, the transmission is much easier to write !
         call MPI_BCAST(eig0_k,nban_dp_rdk*(2*nban_dp_rdk)**formeig ,&
&         MPI_DOUBLE_PRECISION,0,spaceComm,ierr)
       end if
#endif

       if(formeig==0)then
!        The transfer from eig0_k to eig_k uses nban_dp_rdk, which contains
!        the maximal information.
         if(nspinor0==nspinor .or. squeeze==0)then
           eig_k(1:nban_dp_rdk)=eig0_k(1:nban_dp_rdk)
           occ_k(1:nban_dp_rdk)=occ0_k(1:nban_dp_rdk)
         else if(nspinor0==1 .and. nspinor==2)then
           do iband=1,nban_dp_rdk
             eig_k(2*iband  )=eig0_k(iband)
             eig_k(2*iband-1)=eig0_k(iband)
             occ_k(2*iband  )=occ0_k(iband)*0.5_dp
             occ_k(2*iband-1)=occ0_k(iband)*0.5_dp
           end do
         else if(nspinor0==2 .and. nspinor==1)then
           do iband=1,nban_dp_rdk
             eig_k(iband)=eig0_k(2*iband-1)
             occ_k(iband)=occ0_k(2*iband-1)*2.0_dp
           end do
         end if

!        DEBUG
!        write(std_out,*)' wfsinp: me,band_index,ikpt,isppol',me,band_index,ikpt,isppol
!        ENDDEBUG

!        The transfer to eigen uses nban_dp_k, that is bound by the number
!        of bands for this k point.
         ncopy=min(dim_eig_k,(nban_dp_k/nspinor0)*nspinor)
         eigen(1+band_index:ncopy+band_index)=eig_k(1:ncopy)
!        The transfer of occ should be done.

#if defined HAVE_MPI
         if(localrdwf==0 .and. ikpt/=0)then
!          Do not forget : ikpt,isppol were redefined ...
           band_index=band_index_k(ikpt,isppol)
           eigen(1+band_index:(nban_dp_k/nspinor0)*nspinor+band_index) = eig_k(1:(nban_dp_k/nspinor0)*nspinor)
!          The transfer of occ should be done.
         end if
#endif

       else if(formeig==1)then
         call wrtout(std_out,'wfsinp: transfer of first-order eigs not yet coded!',"COLL")
       end if

!      DEBUG
!      write(std_out,*)' wfsinp : me,transferred eig_k',me
!      write(std_out,*)' me,mkmem,nsppol,nsppol0,isppol',me,mkmem,nsppol,nsppol0,isppol
!      write(std_out,*)' me,nkassoc,ikassoc',me,nkassoc,ikassoc
!      ENDDEBUG

       call timab(725,2,tsec)

!      Write to disk if appropriate
!      The coding has to be done ... here, only fragments ...
       call timab(727,3,tsec)

       do isppol_trial=1,nsppol
         if(nsppol==2 .and. nsppol0==2 .and. isppol_trial/=isppol)cycle
         do ikassoc_trial=1,nkassoc

!            DEBUG
!            write(std_out,*)' wfsinp: me, for ikassoc,isppol',&
!            &       me,ikassoc,isppol
!            write(std_out,*)' wfsinp: me, try ikassoc_trial,isppol_trial,nband_k',&
!            &       me,ikassoc_trial,isppol_trial,nband_k
!            ENDDEBUG

!            No conversion is to be done : it will be converted in newkpt
!            If squeeze==0, the block with the ikpt corresponding to ikassoc,
!            and with isppol, contains the wavefunction already
           if( ikassoc_trial/=ikassoc                   .or. &
&           (isppol_trial/=isppol .and. sppoldbl==1 ).or. &
&           squeeze==1                                       )then

             ikpt_trial=indkk0(ikpt0,ikassoc_trial)
             if(sppoldbl==2)then
               if(isppol_trial==1 .and. ikpt_trial>nkpt)cycle
               if(isppol_trial==2 .and. ikpt_trial<=nkpt)cycle
               if(isppol_trial==2)ikpt_trial=ikpt_trial-nkpt
             end if


#if defined HAVE_MPI
             if(ikpt_trial/=0)then
               nband_trial=nband(ikpt_trial+(isppol_trial-1)*nkpt)
               nbd=min(nband_trial,nbd_max);isp=min(isppol_trial,isp_max)
               if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt_trial,1,nbd,isp,me)) ikpt_trial=0
             end if
#endif

             if(ikpt_trial/=0 .and. ikpt_trial<=nkpt_eff)then
               write(message,'(2a,2i5)')ch10,' wfsinp: transfer to ikpt_trial,isppol_trial=',ikpt_trial,isppol_trial
               call wrtout(std_out,message,'PERS')
             end if

             if(ikpt_trial/=0)then
               icg_trial=icg_k(ikpt_trial,isppol_trial)
               band_index_trial=band_index_k(ikpt_trial,isppol_trial)
               nband_trial=nband(ikpt_trial+(isppol_trial-1)*nkpt)
               nban_dp_k=min(nban_dp_rdk,(nband_trial/nspinor)*nspinor0)

               if(squeeze==0)then
!                  GMR: modified to avoid compiler bug
!                  cg(:,1+icg_trial:npw0_k*my_nspinor0*nband_trial+icg_trial)=&
!                  &          cg(:,1+icg:npw0_k*my_nspinor0*nband_trial+icg)
                 do ii=1,npw0_k*my_nspinor0*nband_trial
                   cg(:,ii+icg_trial)=cg(:,ii+icg)
                 end do
!                  GMR
                 if(formeig==0)then
                   eigen(1+band_index_trial:nban_dp_k+band_index_trial)=eig_k(1:nban_dp_k)
!                    occ(1+band_index_trial:nban_dp_k+band_index_trial)=&
!                    &           occ_k(1:nban_dp_k)
                 end if
!                  RF transfer of eigenvalues still to be coded
               else if(squeeze==1)then
                 npw_ktrial=npwarr(ikpt_trial)
                 nband_k=(nban_dp_k/nspinor0)*nspinor
!                  Conversion to be done
                 ceksp=0 ; debug=0 ; icg_disk=0 ; idum=0 ; inplace=0
!                  Note that this routine also convert eig and occ
!                  even if the conversion had already been done


                 call wfconv(ceksp,cg_disk,cg,debug,ecut0,ecut,ecut_eff,&
&                 eig0_k,eig_k,exchn2n3d,formeig,gmet0,gmet,&
&                 icg_disk,icg_trial,ikpt0,ikpt10,ikpt_trial,indkk,&
&                 inplace,isppol_trial,istwfk0,istwfk,&
&                 kg0_k,kg_k,kptns0,kptns,nban_dp_rdk,nband_rdk,&
&                 mcg_disk,mcg,mpi_enreg0,mpi_enreg,mpw0,mpw,&
&                 nban_dp_rdk,nband_trial,ngfft,ngfft,nkpt0,nkpt,&
&                 npw0_k,npw_ktrial,nspinor0,nspinor,nsym,&
&                 occ0_k,occ_k,optorth,randalg,restart,rprimd,&
&                 sppoldbl,symrel,tnons)

!                  DEBUG
!                  write(std_out,*)' wfsinp: ikpt_trial=',ikpt_trial
!                  write(std_out,*)' nband_k,band_index_trial',nband_k,band_index_trial
!                  write(std_out,*)eig0_k(1:nban_dp_k)
!                  write(std_out,*)eig_k(1:nband_k)
!                  ENDDEBUG
                 eigen(1+band_index_trial:nband_k+band_index_trial)=eig_k(1:nband_k)
!                  occ(1+band_index_trial:nband_k+band_index_trial)=&
!                  &           occ_k(1:nband_k)
!                  RF transfer of eigenvalues still to be coded

!                  Endif squeeze==1
               end if

!                DEBUG
!                if(ikpt_trial==2)then
!                write(std_out,*)' wfsinp: iband,ipw,cg for ikpt_trial=2'
!                write(std_out,*)' nband_trial,npw_ktrial=',nband_trial,npw_ktrial
!                do iband=1,nband_trial
!                do ipw=1,npw_ktrial
!                write(std_out,'(2i5,2es16.6)' )&
!                &             iband,ipw,cg(:,ipw+(iband-1)*npw_ktrial+icg_trial)
!                end do
!                end do
!                end if
!                ENDDEBUG

!                End if ikpt_trial/=0
             end if

!              End if ikpt_trial already initialized
           end if

         end do ! ikassoc_trial
       end do ! isppol_trial

       call timab(727,2,tsec)

       ABI_DEALLOCATE(eig_k)
       ABI_DEALLOCATE(eig0_k)
       !if(formeig==0) then
       ABI_DEALLOCATE(occ_k)
       ABI_DEALLOCATE(occ0_k)
       !end if

     end if  ! End the condition of need of this k point
   end do ! End of the k loop
 end do !  End of spin loop

#if defined HAVE_MPI
 call timab(67,1,tsec)

!Still need to skip last k points to read history (MS)
!WARNING : not yet for formeig=1, but is it needed ?
 if(formeig==0)then
   do ikptsp=ikptsp_old+1,nkpt0*nsppol0
     isppol=1 ; if(ikptsp>nkpt0)isppol=2
     ikpt=ikptsp-nkpt0*(isppol-1)
     call WffReadSkipK(formeig,headform0,ikpt,isppol,mpi_enreg,wff1)
   end do
 end if

!Transmit eigenvalues. This routine works in both localrdwf=0 or 1 cases.
 call pareigocc(eigen,formeig,localrdwf,mpi_enreg,mband,nband,nkpt,nsppol,occ,1)


 if(localrdwf==0)then
   ABI_DEALLOCATE(ikpt_me)
   ABI_DEALLOCATE(nband_k_me)
   ABI_DEALLOCATE(ikassoc_me)
   ABI_DEALLOCATE(isppol_me)
 end if

 call timab(67,2,tsec)
#endif

!****************************************************************************

 if(squeeze==1)then
   ABI_DEALLOCATE(kg_k)
   ABI_DEALLOCATE(kg0_k)
 end if


!DEBUG
!if(me==0)then
!write(std_out,*)' wfsinp : cg array='
!icg=0
!do isppol=1,nsppol
!do ikpt=1,1
!nband_k=nband(ikpt+(isppol-1)*nkpt)
!npw_k=npwarr(ikpt)
!do iband=1,nband_k
!write(std_out,*)' new band, icg=',icg
!do ipw=1,npw_k
!write(std_out,'(4i4,2es20.10)' )isppol,ikpt,iband,ipw,cg(:,icg+ipw)
!end do
!icg=icg+npw_k
!end do
!end do
!end do
!end if
!if(ireadwf==1)stop
!write(std_out,*)' wfsinp : eigen array='
!do ikpt=1,nkpt
!do iband=1,mband
!write(std_out,*)'ikpt,iband,eigen',ikpt,iband,eigen(iband+(ikpt-1)*mband)
!end do
!end do
!ENDDEBUG

 ABI_DEALLOCATE(icg_k)
 ABI_DEALLOCATE(band_index_k)

 call timab(720,2,tsec)

end subroutine wfsinp
!!***

!!****f* m_inwffil/initwf
!!
!! NAME
!! initwf
!!
!! FUNCTION
!! Initialization of wavefunctions.
!! If formeig==1, and partially filled case, I am not sure that the eig_k are initialized properly ...
!! formeig option (format of the eigenvalues and eigenvector) :
!!   0 => ground-state format (initialisation of
!!        eigenvectors with random numbers, vector of eigenvalues)
!!   1 => respfn format (initialisation of
!!        eigenvectors with 0 s, hermitian matrix of eigenvalues)
!!
!! INPUTS
!! formeig=see above
!! headform=header format (might be needed to read the block of wfs)
!! icg=shift to be given to the location of the data in the array cg
!! ikpt= number of the k point of which the wf is initialised
!! spin=spin index
!! mcg=dimension of the cg array
!! mpi_enreg=information about MPI parallelization
!! nband_k=number of bands at this particular k point
!! nkpt=number of k points
!! npw=number of plane waves
!! nspinor=number of spinorial components of the wavefunctions (on current proc)
!! wff1=structure info for file containing wavefunctions (when needed)
!!
!! OUTPUT
!! cg(2,mcg)=complex wf array
!!  if ground state format (formeig=0):
!!    eig_k(nband_k)=list of eigenvalues (input or init to large number), hartree
!!  if respfn format (formeig=1):
!!    eig_k(2*nband_k*nband_k)= matrix of eigenvalues (input or init to large number), hartree
!!
!! SIDE EFFECTS
!! Input/output:
!! occ_k(nband_k)=list of occupations (input or left to their initial value)
!! ikptsp_old=number of the previous spin-k point, or 0 if first call of present file
!!
!! PARENTS
!!      wfsinp
!!
!! CHILDREN
!!      rwwf,timab,wffreadskipk,wfk_close,wfk_open_read,wfk_read_band_block
!!      wrtout
!!
!! SOURCE

subroutine initwf(cg,eig_k,formeig,headform,icg,ikpt,ikptsp_old,&
&  spin,mcg,mpi_enreg,nband_k,nkpt,npw,nspinor,occ_k,wff1)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: formeig,headform,icg,ikpt,spin,mcg,nband_k,nkpt,npw,nspinor
 integer,intent(inout) :: ikptsp_old
 type(MPI_type),intent(in) :: mpi_enreg
 type(wffile_type),intent(inout) :: wff1
!arrays
 real(dp),intent(inout) :: occ_k(nband_k)
 real(dp),intent(inout) :: cg(2,mcg),eig_k((2*nband_k)**formeig*nband_k) !vz_i

!Local variables-------------------------------
!scalars
 integer,parameter :: nkpt_max=50
 integer :: ikpt0,nband_disk,tim_rwwf
 character(len=500) :: msg
!arrays
 integer,allocatable :: kg_dum(:,:)
 real(dp) :: tsec(2)
#if 0
 integer :: iomode,comm,funt,ierr
 type(wfk_t) :: Wfk
#endif

! *************************************************************************

!DEBUG
!write(std_out,*)' initwf : enter, ikptsp_old,ikpt,spin,nkpt= ',ikptsp_old,ikpt,spin,nkpt
!stop
!ENDDEBUG

#if 0
 MSG_WARNING("Entering new IO section")
!call WffClose(wff1,ierr)
 comm   = MPI_enreg%comm_cell
 iomode = iomode_from_fname(wff1%fname)
 call wfk_open_read(Wfk,wff1%fname,formeig,iomode,get_unit(),comm)
 call wfk_read_band_block(Wfk,(/1,nband_k/),ikpt,spin,xmpio_at,cg_k=cg(1:,icg:),eig_k=eig_k,occ_k=occ_k)
 call wfk_close(Wfk)
!call clsopn(wff1)
 RETURN
#endif

 call timab(770,1,tsec)
 call timab(771,1,tsec)

 ABI_ALLOCATE(kg_dum,(3,0))

!Skip wavefunctions for k-points not treated by this proc.
!(from ikptsp_old+1 to ikpt+(spin-1)*nkpt-1)
 if (ikptsp_old<ikpt+(spin-1)*nkpt-1) then
   do ikpt0=ikptsp_old+1,ikpt+(spin-1)*nkpt-1
     call WffReadSkipK(formeig,headform,ikpt0,spin,mpi_enreg,wff1)
   end do
 end if

!DEBUG
!write(std_out,*)' initwf : before rwwf'
!write(std_out,*)' formeig,icg,ikpt,spin=',formeig,icg,ikpt,spin
!write(std_out,*)' nband_k,nband_disk,npw,nspinor=',nband_k,nband_disk,npw,nspinor
!write(std_out,*)' unwff1=',unwff1
!stop
!ENDDEBUG

 if(mpi_enreg%paralbd==0)tim_rwwf=2
 if(mpi_enreg%paralbd==1)tim_rwwf=20

 call timab(771,2,tsec)

 call rwwf(cg,eig_k,formeig,headform,icg,ikpt,spin,kg_dum,nband_k,mcg,mpi_enreg,nband_k,nband_disk,&
& npw,nspinor,occ_k,1,0,tim_rwwf,wff1)

 call timab(772,1,tsec)

 if (ikpt<=nkpt_max) then
   write(msg,'(3(a,i0))')' initwf: disk file gives npw= ',npw,' nband= ',nband_disk,' for kpt number= ',ikpt
   call wrtout(std_out,msg,'PERS')
 else if (ikpt==nkpt_max+1) then
   call wrtout(std_out,' initwf: the number of similar message is sufficient... stop printing them','PERS')
 end if

!Check the number of bands on disk file against desired number. These are not required to agree)
 if (nband_disk /= nband_k) then
   write(msg,'(2(a,i0),3a,i0,3a)')&
   'For kpt number ',ikpt,' disk file has ',nband_disk,' bands',ch10,&
   'but input file gave nband= ',nband_k,'.',ch10,&
   'This is not fatal. Bands are skipped or filled with random numbers.'
   MSG_COMMENT(msg)
 end if

 if (ikpt<=nkpt_max) then
   write(msg,'(a,i0,a)')' initwf: ',nband_disk,' bands have been initialized from disk'
   call wrtout(std_out,msg,'PERS')
 end if

 ikptsp_old=ikpt+(spin-1)*nkpt

 ABI_DEALLOCATE(kg_dum)

 call timab(772,2,tsec)
 call timab(770,2,tsec)

end subroutine initwf
!!***

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
!!  mpi_enreg1=information about MPI parallelization, for the input wf file
!!  mpi_enreg2=information about MPI parallelization, for the output wf file
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

subroutine newkpt(ceksp2,cg,debug,ecut1,ecut2,ecut2_eff,eigen,exchn2n3d,fill,&
&                  formeig,gmet1,gmet2,headform1,indkk,iout,ireadwf,&
&                  istwfk1,istwfk2,kg2,kptns1,kptns2,mband2,mcg,mkmem1,mkmem2,&
&                  mpi_enreg1,mpi_enreg2,mpw1,mpw2,my_nkpt2,nband1,nband2,&
&                  ngfft1,ngfft2,nkpt1,nkpt2,npwarr1,npwarr2,nspinor1,nspinor2,&
&                  nsppol1,nsppol2,nsym,occ,optorth,prtvol,randalg,restart,rprimd,&
&                  sppoldbl,symrel,tnons,unkg2,wffinp,wffout)

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
   write(message, '(5a,i2,a,i2,2a,i2,a,i2,4a)' ) &
&   'The wavefunction translator is (still) unable to interchange',ch10,&
&   'spin-polarized wfs and spinor wfs. However,',ch10,&
&   'the input  variables are nsppol1=',nsppol1,', and nspinor1=',nspinor1,ch10,&
&   'the output variables are nsppol2=',nsppol2,', and nspinor2=',nspinor2,ch10,&
&   'Action: use a non-spin-polarized wf to start a spinor wf,',ch10,&
&   '        and a non-spinor wf to start a spin-polarized wf.'
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
&     'The dimension mcg= ',mcg,', should be larger than',ch10,&
&     'the product of mband2= ',mband2,', mpw2= ',mpw2,', and nspinor2= ',my_nspinor2
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
         call wrtout(std_out, '- newkpt: prtvol=0 or 1, do not print more k-points.', 'PERS')
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
         ' newkpt: about to call randac',ch10,&
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
&         'When mkmem1/=0, one must have nbd1=nbd1_rd, while',ch10,&
&         'nbd1 = ',nbd1,', and nbd1_rd = ',nbd1_rd
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
&         ' newkpt: about to call rwwf with ikpt1=',ikpt1,ch10,&
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

!!****f* ABINIT/wfconv
!! NAME
!! wfconv
!!
!! FUNCTION
!! This subroutine treats the wavefunctions for one k point,
!! and converts them to other parameters.
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
!!  mpi_enreg1=information about MPI parallelization for set 1
!!  mpi_enreg2=information about MPI parallelization for set 2
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

subroutine wfconv(ceksp2,cg1,cg2,debug,ecut1,ecut2,ecut2_eff,&
& eig_k1,eig_k2,exchn2n3d,formeig,gmet1,gmet2,icg1,icg2,&
& ikpt1,ikpt10,ikpt2,indkk,inplace,isppol2,istwfk1,istwfk2,&
& kg1,kg2,kptns1,kptns2,mband1,mband2,mcg1,mcg2,mpi_enreg1,mpi_enreg2,&
& mpw1,mpw2,nbd1,nbd2,ngfft1,ngfft2,nkpt1,nkpt2,npw1,npw2,nspinor1,nspinor2,&
& nsym,occ_k1,occ_k2,optorth,randalg,restart,rprimd2,sppoldbl,symrel,tnons)

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
 !write(std_out,*)' wfconv : isym=',isym
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
       write(message,'(a,i6,a,i7,a,i4)')' wfconv:',nbremn,' bands set=0 with npw=',npw2,', for ikpt=',ikpt2
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

end module m_inwffil
!!***
