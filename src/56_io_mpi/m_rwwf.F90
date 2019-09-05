!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_rwwf
!! NAME
!!  m_rwwf
!!
!! FUNCTION
!!   Read/Write wavefunctions.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA,XG,GMR,MVer,MB,MT)
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

module m_rwwf

 use defs_basis
 use m_errors
 use m_wffile
 use m_abicore
 use m_xmpi
#if defined HAVE_MPI2
 use mpi
#endif
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_abitypes, only : mpi_type
 use m_time,   only : timab

 implicit none

 private

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!!***

 public :: rwwf
 public :: WffReadSkipK
!!***

contains
!!***

!!****f* m_rwwf/rwwf
!! NAME
!! rwwf
!!
!! FUNCTION
!!  This subroutine reads (different options) or write (option=2) the block of records
!!  related to one k point, and one spin-polarization, that
!!  contains the wavefunctions (as well as the eigenvalues and occupations).
!!  If called with option -1, the records will be skipped.
!!  If called with option -2, only the wavefunctions are read.
!!  The disk file unitwf should have been prepared
!!  outside of this routine, in order to read or write the correct records.
!!
!! INPUTS
!!  formeig=format of the eigenvalues
!!     0 => vector of eigenvalues
!!     1 => hermitian matrix of eigenvalues
!!  headform=format of the header of the wf file, also governing the k block format
!!    in case headform=0, use the default (current) format and headform
!!  icg=shift to be given to the location of the cg array
!!  ikpt=index of current k point (only needed for error message)
!!  isppol=spin polarization currently treated (only needed for error message)
!!  mband=maximum number of bands (dimension of cg, eigen and occ)
!!  mcg=dimention of cg
!!  nband=number of bands actually in cg, eigen and occ
!!   (if writing mode : must be larger or equal to nband_disk, only nband_disk bands are written ;
!!    if reading mode : can be equal, larger or smaller than nband_disk, but
!!     cg, eigen and occ will not be completely filled if nband>nband_disk)
!!  nband_disk=number of bands on the disk file
!!  npw=number of plane waves
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  option= 2 for writing cg, eigen and occ,
!!          1 for reading cg and eigen,
!!         -1 for reading/skipping,
!!         -2 for reading cg only
!!          3 for reading the eigenvalues only
!!          4 for writing a file containing only eigenvalues and occupations (need to be read with option 4)
!!          5 for writing a file containing only eigenvalues and occupations, that can however be read with option 3)
!!         -4 for reading a file written with 4
!!                 (different from 3 which reads a normal option 2 file)
!!  optkg= if 1 , read or write kg_k ; if 0, do not care about kg_k
!!  tim_rwwf=timing code of the calling routine (set to 0 if not attributed)
!!  wff=struct info for wavefunction
!!   | unitwf=unit number for wavefunction
!!
!! SIDE EFFECTS
!!  cg(2,npw*nspinor*mband)=planewave coefficients of wavefunctions,
!!    input if option=2; output if option=1 or -2
!!  eigen((2*mband)**formeig *mband)=array for holding eigenvalues (hartree)
!!    input if option=2 or 4 or 5; output if option=1
!!  kg_k(3,optkg*npw)=k+g data  (only if optkg==1)
!!    input if option=2; output if option=1 or -2
!!  nband_disk=number of bands on disk
!!    input if option=2 or 4 or 5; output in the other cases
!!  occ(mband)=array for holding eigenvalues (hartree)
!!    input if option=2 or 4 or 5; output if option=1
!!    no meaning if frmeig/=0
!!
!! NOTES
!!  WARNING : occ is not read in the present status of this routine
!!  WARNING : skipping k-blocks is also done in the randac subroutine
!!  WARNING : reading the two first records is also done in the rdnpw routine
!!  WARNING : writing the two first records is also done in the dfpt_vtowfk routine
!!
!! TODO!     if (mpi_enreg%flag_ind_kg_mpi_to_seq==1 ) then
!!  Some arguments are contained in the wff datastructure, and should be eliminated.
!!  option 3 should be called -3 (reading -> negative option) and others (-1,1) re-shuffled.
!!
!! PARENTS
!!      WffReadEigK,WffReadSkipK,initwf,m_iowf,m_wfk,newkpt,uderiv
!!
!! CHILDREN
!!      mpi_bcast,wffreadwrite_mpio,wffwritenpwrec,xderivewrecend
!!      xderivewrecinit,xderivewrite,xmpi_sum
!!
!! SOURCE

subroutine rwwf(cg,eigen,formeig,headform,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,&
&               nband,nband_disk,npw,nspinor,occ,option,optkg,tim_rwwf,wff)

!Arguments ------------------------------------
 integer,intent(in) :: formeig,headform,icg,ikpt,isppol,mband,mcg,nband,npw
 integer,intent(inout) :: nband_disk
 integer,intent(in) :: nspinor,option,optkg,tim_rwwf
 integer,intent(inout),target :: kg_k(3,optkg*npw)
 real(dp),intent(inout),target :: cg(2,mcg),eigen((2*mband)**formeig*mband),occ(mband)
 type(wffile_type),intent(inout) :: wff
 type(MPI_type), intent(in) :: mpi_enreg

!Local variables-------------------------------
!scalars
 character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)

! *************************************************************************

 call timab(270+tim_rwwf,1,tsec)

!Might check that icg+npw*nband*nspinor is smaller than mcg

!Check that nband is smaller than mband, if one will not skip the records.
 if (nband>mband .and. option/=-1)then
   write(msg,'(a,i0,a,i0,a)')' One should have nband<=mband. However, nband=',nband,', and mband=',mband,'.'
   MSG_BUG(msg)
 end if

!Check that formeig is 0 or 1.
 if ( ALL(formeig/=(/0,1/)) ) then
   write(msg,'(a,i0,a)')' The argument formeig should be 0 or 1. However, formeig=',formeig,'.'
   MSG_BUG(msg)
 end if

!Check the value of option
 if ( ALL( option /= (/1,2,3,4,5,-1,-2,-4/) )) then
   write(msg,'(a,i0,a)')' The argument option should be between 1 and 5, or -1, -2, -4. However, option=',option,'.'
   MSG_BUG(msg)
 end if

 if (option/=2.and.option/=4 .and. option/=5 ) then ! read
   call readwf(cg,eigen,formeig,headform,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,&
&   nband,nband_disk,npw,nspinor,occ,option,optkg,wff)
 else                                ! write
   call writewf(cg,eigen,formeig,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,&
&   nband,nband_disk,npw,nspinor,occ,option,optkg,wff)
 end if

 call timab(270+tim_rwwf,2,tsec)

end subroutine rwwf
!!***

!!****f* m_rwwf/WffReadSkipK
!! NAME
!! WffReadSkipK
!!
!! FUNCTION
!!  (Wavefunction file, read action : skip one k-point blok)
!!  This subroutine skips the block of records
!!  related to one k point, and one spin-polarization, that
!!  contains the wavefunctions as well as the eigenvalues and occupations,
!!  in a wavefunction file that has been already initialized.
!!
!! INPUTS
!!  formeig=format of the eigenvalues
!!   0 => vector of eigenvalues (for Ground-State files)
!!   1 => hermitian matrix of eigenvalues (for Response-Function files)
!!  headform=format of the header of the wf file, also governing the k block format
!!   in case headform=0, use the default (current) format and headform
!!  ikpt=index of current k point (only needed for error message)
!!  isppol=spin polarization currently treated (only needed for error message)
!!  mpi_enreg=information about MPI parallelization
!!  wff=structured info for wavefunction file
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      initwf,newkpt,wfsinp
!!
!! CHILDREN
!!      rwwf
!!
!! SOURCE

subroutine WffReadSkipK(formeig,headform,ikpt,isppol,mpi_enreg,wff)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: formeig,headform,ikpt,isppol
 type(MPI_type),intent(in) :: mpi_enreg
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
!scalars
 integer :: icg,mband,mcg,nband,nband_disk,option,optkg,tim_rwwf
 integer,parameter :: nspinor1=1,npw1=1
!arrays
 integer,allocatable :: kg_dum(:,:)
 real(dp) :: cg_dum(2,1),occ_dum(1)
 real(dp),allocatable :: eig_dum(:)

! *************************************************************************

 ! No need to skip if netcdf
 if (wff%iomode == IO_MODE_ETSF) return

 option=-1
 tim_rwwf=0 ; mcg=1 ; mband=1 ; icg=0 ; optkg=0 ; nband=0
 ABI_ALLOCATE(eig_dum,(2**formeig))
 ABI_ALLOCATE(kg_dum,(3,optkg*npw1))

 call rwwf(cg_dum,eig_dum,formeig,headform,icg,ikpt,isppol,kg_dum,mband,mcg,mpi_enreg,nband,&
& nband_disk,npw1,nspinor1,occ_dum,option,optkg,tim_rwwf,wff)

 ABI_DEALLOCATE(eig_dum)
 ABI_DEALLOCATE(kg_dum)

end subroutine WffReadSkipK
!!***

! -------------------------------------------------------------------------------------------------

!!****f* m_rwwf/readwf
!! NAME
!! readwf
!!
!! FUNCTION
!!  This subroutine reads the block of records related to one k point, and one spin-polarization, that
!!  contains the wavefunctions (as well as the eigenvalues).
!!  The disk file unitwf should have been prepared outside of this routine, in order to read the correct records.
!!
!! INPUTS
!!  formeig=format of the eigenvalues
!!    0 => vector of eigenvalues
!!    1 => hermitian matrix of eigenvalues
!!  icg=shift to be given to the location of the cg array
!!  ikpt=index of current k point (only needed for error message)
!!  isppol=spin polarization currently treated (only needed for error message)
!!  kg_k(3,optkg*npw)=k+g data  (only if optkg==1)
!!  mband=maximum number of bands (dimension of cg, eigen and occ)
!!  mcg=dimension of cg
!!  nband=number of bands actually in cg, eigen and occ
!!  nband_disk=number of bands on the disk file
!!  npw=number of plane waves (on current proc)
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  option= 1 for reading cg, eigen and occ,
!!         -1 for reading/skipping,
!!         -2 for reading cg only
!!          3 for reading the eigenvalues only
!!         -4 for reading a file written with 4
!!  optkg= if 1 , read or write kg_k ; if 0, do not care about kg_k
!!  wff=struct info for wavefunction
!!
!! SIDE EFFECTS
!!  Current kpt and spin updated
!!  cg(2,npw*nspinor*mband)=planewave coefficients of wavefunctions,
!!  eigen((2*mband)**formeig *mband)=array for holding eigenvalues (hartree)
!!  occ(mband)=array for holding electronic occupations
!!
!! NOTES
!!  WARNING : occ is not read in the present status of this routine
!!  WARNING : skipping k-blocks is also done in the randac subroutine
!!  WARNING : reading the two first records is also done in the rdnpw routine
!!
!! PARENTS
!!      rwwf
!!
!! CHILDREN
!!      mpi_bcast,wffreadwrite_mpio,wffwritenpwrec,xderivewrecend
!!      xderivewrecinit,xderivewrite,xmpi_sum
!!
!! SOURCE

subroutine readwf(cg,eigen,formeig,headform,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,&
&                 nband,nband_disk,npw,nspinor,occ,option,optkg,wff)

!Arguments ------------------------------------
 integer,intent(in) :: formeig,headform,icg,ikpt,isppol,mband,mcg,nband,npw
 integer,intent(in) :: nspinor,option,optkg
 integer,intent(inout) :: nband_disk
 integer,intent(inout),target :: kg_k(3,optkg*npw)
 real(dp),intent(inout),target :: cg(2,mcg),eigen((2*mband)**formeig*mband),occ(mband)
 type(MPI_type),intent(in) :: mpi_enreg
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
 integer :: iband,indxx,ios,ipw,ispinor_index,nband1,ncid_hdr
 integer :: npw1,npwso,npwso1,npwsotot,npwtot,nrec,nspinor1,nspinortot,unitwf,use_f90
 integer :: band1,band2,ierr
 character(len=500) :: msg
 character(len=fnlen) :: fname
 integer :: ikpt_this_proc,ispinor
 integer,allocatable :: ind_cg_mpi_to_seq(:)
#ifdef HAVE_NETCDF
 integer :: kg_varid,eig_varid,occ_varid,cg_varid,h1_varid,mband_varid,ncerr,ii,mband_file
 integer,allocatable :: gall(:,:)
 real(dp),allocatable :: cg_all(:,:,:),h1mat(:,:,:)
#endif

! *********************************************************************

 !write(std_out,*)"rwwf with ikpt, isppol, option, etsf",ikpt, isppol, option, wff%iomode == IO_MODE_ETSF

!Check the options
 if ( ALL(option /= [1,3,-1,-2,-4])) then
   write(msg,'(a,i0,a)')'The argument option should be -4, -2, -1, 1 or 3.  However, option=',option,'.'
   MSG_BUG(msg)
 end if

 npwtot=npw; npwso=npw*nspinor
 npwsotot=npwso
 nspinortot=min(2,(1+mpi_enreg%paral_spinor)*nspinor)

 unitwf=wff%unwff
 ncid_hdr=unitwf

 use_f90=0
 if (wff%iomode==IO_MODE_FORTRAN.or.(wff%iomode==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) use_f90=1

 if (option==1.or.option==-2) then

   ! Compute mapping my_gtable --> sequential gtable.
   if (any(wff%iomode==[IO_MODE_MPI, IO_MODE_ETSF]).and.nband>0) then
     call xmpi_sum(npwsotot,wff%spaceComm_mpiio,ios)
     npwtot=npwsotot/nspinortot
     ABI_ALLOCATE(ind_cg_mpi_to_seq,(npwso))
     if (allocated(mpi_enreg%my_kgtab)) then
       ikpt_this_proc=mpi_enreg%my_kpttab(ikpt)
       if ( ikpt_this_proc <= 0  ) then
         MSG_BUG("rwwf: ikpt_this_proc <= 0")
       end if
       do ispinor=1,nspinor
         ispinor_index=ispinor
         if (mpi_enreg%nproc_spinor>1) ispinor_index=mpi_enreg%me_spinor + 1
         ind_cg_mpi_to_seq(1+npw*(ispinor-1):npw*ispinor)=npwtot*(ispinor_index-1) &
&         + mpi_enreg%my_kgtab(1:npw,ikpt_this_proc)
       end do
     else
       ind_cg_mpi_to_seq(1:npwso) = (/(ipw,ipw=1,npwso)/)
     end if
   end if
 end if

!---------------------------------------------------------------------------
!Read the first record: npw, nspinor, nband_disk
!---------------------------------------------------------------------------

 if (headform>=40.or.headform==0) then ! headform==0 refers to the current headform
   call WffReadNpwRec(ios,ikpt,isppol,nband_disk,npw1,nspinor1,wff)
   npwso1=npw1*nspinor1
   if(ios/=0)then
     inquire (unit=unitwf, NAME=fname)
     write(msg,'(3a,i4,2a,i4,4a)') &
&     'Reading option of rwwf. Trying to read',ch10,&
&     'the (npw,nspinor,nband) record of a wf file, unit=',unitwf,ch10,&
&     'gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&     'Action: check your input wf file:',trim(fname)
     MSG_ERROR(msg)
   end if

 else
   ! Old format
   if(use_f90==1)then
     read (unitwf,iostat=ios) npwso1,nband_disk
   else if(wff%iomode==IO_MODE_MPI)then
     call xderiveRRecInit(wff,ios)
     call xderiveRead(wff,npwso1,ios)
     call xderiveRead(wff,nband_disk,ios)
     call xderiveRRecEnd(wff,ios)
   end if
   if(ios/=0)then
     inquire (unit=unitwf, NAME=fname)
     write(msg,'(3a,i4,2a,i4,4a)') &
&     'Reading option of rwwf. Trying to read',ch10,&
&     'the (npw,nband) record of a wf file with headform <40 , unit=',unitwf,ch10,&
&     'gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&     'Action: check your input wf file:',trim(fname)
     MSG_ERROR(msg)
   end if
 end if ! headform

 if (option==1.or.option==-2) then !  Will read the wavefunction and/or kg data, so check npw and nspinor

   if (headform>=40.or.headform==0) then ! New format. headform==0 refers to the current headform
     if (npwtot/=npw1) then
       write(msg,'(3a,i0,a,i0,a)') &
&       'Reading option of rwwf. One should have npwtot=npw1',ch10,&
&       'However, npwtot= ',npwtot,', and npw1= ',npw1,'.'
       MSG_BUG(msg)
     end if
     if(nspinortot/=nspinor1)then
       write(msg,'(3a,i0,a,i0,a)') &
&       'Reading option of rwwf. One should have nspinor=nspinor1',ch10,&
&       'However, nspinortot= ',nspinortot,', and nspinor1= ',nspinor1,'.'
       MSG_BUG(msg)
     end if
   else ! Treat the Old format.
     if(npwsotot/=npwso1)then
       write(msg,'(3a,i0,a,i0,a)') &
&       'Reading option of rwwf. One should have npwso=npwso1',ch10,&
&       'However, npwsotot= ',npwsotot,', and npwso1= ',npwso1,'.'
       MSG_BUG(msg)
     end if
   end if ! headform
!
 end if ! option==1.or.option==2

!---------------------------------------------------------------------------
!Read the second record: (k+G) vectors
!---------------------------------------------------------------------------

 if (headform>=40.or.headform==0) then ! headform==0 refers to the current headform

   if ((option==1.or.option==-2.or.option==3).and.optkg/=0 )then

     if(use_f90==1)then
       read(unitwf,iostat=ios) kg_k(1:3,1:npw)

     else if(wff%iomode==IO_MODE_MPI)then
       call xderiveRRecInit(wff,ios)
       if (allocated(mpi_enreg%my_kgtab)) then
         ikpt_this_proc=mpi_enreg%my_kpttab(ikpt)
         call xderiveRead(wff,kg_k(1:3,1:npw),3,npw,wff%spaceComm_mpiio,&
&         mpi_enreg%my_kgtab(1:npw,ikpt_this_proc),ios)
       else
!        MG The call below uses MPI_SCAN but here we want to read the full set of G.
         call xderiveRead(wff,kg_k(1:3,1:npw),3,npw,wff%spaceComm_mpiio,ios)
         call xmpi_read_int2d(wff,kg_k(1:3,1:npw),wff%spaceComm_mpiio,xmpio_collective,ios)
       end if
       call xderiveRRecEnd(wff,ios)

#ifdef HAVE_NETCDF
     else if (wff%iomode == IO_MODE_ETSF) then
       ! Read reduced_coordinates_of_plane_waves for this k point (npw1 is npw_disk).
       ! TODO: spinor parallelism
       NCF_CHECK(nf90_inq_varid(wff%unwff, "reduced_coordinates_of_plane_waves", kg_varid))
       if (npw == npw1) then
         ncerr = nf90_get_var(wff%unwff, kg_varid, kg_k, start=[1,1,ikpt], count=[3,npw1,1])
         NCF_CHECK_MSG(ncerr, "getting kg_k")
       else
         write(std_out,*)"ETSF Reading distributed kg_k"
         ABI_MALLOC(gall, (3, npw1))
         ncerr = nf90_get_var(wff%unwff, kg_varid, gall, start=[1,1,ikpt], count=[3,npw1,1])
         NCF_CHECK_MSG(ncerr, "getting kg_k")
         do ipw=1,npw
           kg_k(:,ipw) = gall(:,ind_cg_mpi_to_seq(ipw))
         end do
         ABI_FREE(gall)
       end if
#endif
     end if

   else ! option
     call WffReadSkipRec(ios,1,wff) ! Skip the record
   end if

   if(ios/=0)then
     inquire (unit=unitwf, NAME=fname)
     write(msg,'(3a,i4,2a,i4,4a)')  &
&     'Reading option of rwwf. Trying to read',ch10,&
&     'the k+g record of a wf file, unit=',unitwf,ch10,&
&     'gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&     'Action: check your input wf file:',trim(fname)
     MSG_ERROR(msg)
   end if
 end if ! headform

!---------------------------------------------------------------------------
!Read the third record: eigenvalues
!---------------------------------------------------------------------------
!The reading of occ should be enabled, BUT taking into account
!of headform of the disk file : occ was NOT present in the disk files with headform=22

 nband1 = min(nband,nband_disk)

!===== Case formeig=0: read eigenvalues =====
 if (formeig==0) then

   if (option==1.or.option==3.or.option==-4) then
     if(use_f90==1)then
       read (unitwf,iostat=ios) eigen(1:nband1)
     else if(wff%iomode==IO_MODE_MPI)then
       call xderiveRRecInit(wff,ios)
       call xderiveRead(wff,eigen,nband1,xmpi_comm_self,ios)
       call xderiveRRecEnd(wff,ios)
     end if
   else
     call WffReadSkipRec(ios,1,wff)
   end if ! option

   if(ios/=0)then
     inquire (unit=unitwf, NAME=fname)
     write(msg,'(3a,i4,2a,i4,4a)') &
&     'Reading option of rwwf. Trying to read',ch10,&
&     'an eigenvalue record of a wf file, unit=',unitwf,ch10,&
&     'gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&     'Action: check your input wf file.',trim(fname)
     MSG_ERROR(msg)
   end if

#ifdef HAVE_NETCDF
   if (wff%iomode == IO_MODE_ETSF) then
     ! get eigenvalues and occupations
     NCF_CHECK(nf90_inq_varid(wff%unwff, "eigenvalues", eig_varid))
     ncerr = nf90_get_var(wff%unwff, eig_varid, eigen, start=[1,ikpt,isppol], count=[nband1,1,1])
     NCF_CHECK_MSG(ncerr, "getting eig_k")

     NCF_CHECK(nf90_inq_varid(wff%unwff, "occupations", occ_varid))
     ncerr = nf90_get_var(wff%unwff, occ_varid, occ, start=[1,ikpt,isppol], count=[nband1,1,1])
     NCF_CHECK_MSG(ncerr, "getting occ_k")
   end if
#endif

!  ===== Case formeig=1: read matrix of eigenvalues =====
!  Will be written later (together with wave-functions)
 else if(formeig==1)then
 end if ! formeig

!---------------------------------------------------------------------------
!Read the wave-function coefficients
!---------------------------------------------------------------------------

!Select bands
 nband1=min(nband,nband_disk)
 if(nband1>0.and.option/=-1)then

!  ===== Case formeig=0: read only wave-functions =====
   if (formeig==0) then

     if (option==1.or.option==-2) then

       if (use_f90==1) then
         do iband=1,nband1
           ipw=(iband-1)*npwso+icg
           read(unitwf,iostat=ios) cg(1:2,ipw+1:ipw+npwso)
           if (ios/=0) exit
         end do
       else if (wff%iomode==IO_MODE_MPI) then
         call WffReadWrite_mpio(wff,1,cg,mcg,icg,nband1,npwso,npwsotot,ind_cg_mpi_to_seq,ios)
       end if

     else if (option/=-4) then
       do iband=1,nband1
         call WffReadSkipRec(ios,1,wff) ! Skip the record
       end do
     end if ! option

     if(ios/=0)then
       inquire (unit=unitwf, NAME=fname)
       write(msg,'(3a,i4,2a,i4,4a)') &
&       'Reading option of rwwf. Trying to read',ch10,&
&       'a RF wf record of a wf file, unit=',unitwf,ch10,&
&       'gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&       'Action: check your input wf file.',trim(fname)
       MSG_ERROR(msg)
     end if

#ifdef HAVE_NETCDF
     if (wff%iomode == IO_MODE_ETSF) then
       ! The coefficients_of_wavefunctions on file have shape [cplex, mpw, nspinor, mband, nkpt, nsppol]
       NCF_CHECK(nf90_inq_varid(wff%unwff, "coefficients_of_wavefunctions", cg_varid))
       if (npw == npw1) then
         ncerr = nf90_get_var(wff%unwff, cg_varid, cg(:, icg+1:), start=[1,1,1,1,ikpt,isppol], &
         count=[2,npw,nspinor,nband1,1,1])
         NCF_CHECK_MSG(ncerr, "getting cg_k")
       else
         write(std_out,*)"ETSF Reading distributed cg"
         ABI_MALLOC_OR_DIE(cg_all, (2, npw1*nspinor, nband1), ierr)
         ncerr = nf90_get_var(wff%unwff, cg_varid, cg_all, start=[1,1,1,1,ikpt,isppol], &
         count=[2,npw1,nspinor,nband1,1,1])
         NCF_CHECK_MSG(ncerr, "getting cg_k")
         ii = icg
         do iband=1,nband1
           do ipw=1,npw
             ii = ii + 1
             cg(:,ii) = cg_all(:,ind_cg_mpi_to_seq(ipw),iband)
           end do
         end do
         ABI_FREE(cg_all)
       end if
     end if
#endif

!    ===== Case formeig=1: read eigenvalues, occupations and wave-functions =====
   else if(formeig==1)then
     ! write(std_out,*)"nband1",nband1
     ! ABI_CHECK(nband1==nband_disk,"nband != nband_disk")

     if (wff%iomode /= IO_MODE_ETSF) then
       indxx=0
       do iband=1,nband1

         if (option==1.or.option==3.or.option==-4) then
           if(use_f90==1)then
             read (unitwf,iostat=ios) eigen(1+indxx:2*nband1+indxx)
           else if(wff%iomode==IO_MODE_MPI)then
!            Should use an access with a "view"
             call xderiveRRecInit(wff,ios)
             call xderiveRead(wff,eigen(1+indxx:2*nband1+indxx),2*nband1,xmpi_comm_self,ios)
             call xderiveRRecEnd(wff,ios)
           end if
           indxx=indxx+2*nband1
         else
           call WffReadSkipRec(ios,1,wff) ! Skip the record
         end if

         if(ios/=0)then
           inquire (unit=unitwf, NAME=fname)
           write(msg,'(3a,i4,2a,i4,4a)') &
&           'Reading option of rwwf. Trying to read',ch10,&
&           'a RF eigenvalue record of a wf file, unit=',unitwf,ch10,&
&           'gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&           'Action: check your input wf file.',trim(fname)
           MSG_ERROR(msg)
         end if

         if(option==1.or.option==-2)then
           ipw=(iband-1)*npwso+icg
           if(use_f90==1)then
             ipw=(iband-1)*npwso+icg
             read(unitwf,iostat=ios) cg(1:2,ipw+1:ipw+npwso)
           else if(wff%iomode==IO_MODE_MPI)then
!            Should use an access with a "view"
             call xderiveRRecInit(wff,ios)
             call xderiveRead(wff,cg(1:2,ipw+1:ipw+npwso),2,npwso,wff%spaceComm_mpiio,ind_cg_mpi_to_seq,ios)
             call xderiveRRecEnd(wff,ios)
           end if
         else if (option/=-4) then
           call WffReadSkipRec(ios,1,wff) ! Skip the record
         end if ! option

         if(ios/=0)then
           inquire (unit=unitwf, NAME=fname)
           write(msg,'(3a,i4,2a,i4,4a)') &
&           'Reading option of rwwf. Trying to read',ch10,&
&           'a RF wf record of a wf file, unit=',unitwf,ch10,&
&           'gave iostat=',ios,'. Your file is likely not correct.',ch10,&
&           'Action: check your input wf file.',trim(fname)
           MSG_ERROR(msg)
         end if
       end do ! iband

     else
#ifdef HAVE_NETCDF
        ! ETSF-IO
       if (any(option == [1, 3, -4])) then
          ! Read eigen. Remember that the matrix on file has shape [2, mband, mband, nkpt, nspin]
          ! whereas the eigen array used in Abinit is packed. Read the full matrix first, then pack data.
         NCF_CHECK(nf90_inq_dimid(wff%unwff, "max_number_of_states", mband_varid))
         NCF_CHECK(nf90_inquire_dimension(wff%unwff, mband_varid, len=mband_file))
         h1_varid = nctk_idname(wff%unwff, "h1_matrix_elements")

         ABI_MALLOC(h1mat, (2, mband_file, mband_file))
         ncerr = nf90_get_var(wff%unwff, h1_varid, h1mat, start=[1,1,1,ikpt,isppol])
         NCF_CHECK_MSG(ncerr, "getting h1_matrix_elements")

         indxx=1
         do band2=1,nband1
           do band1=1,nband1
             eigen(indxx:indxx+1) = h1mat(:,band1,band2)
             indxx = indxx + 2
           end do
         end do
         ABI_FREE(h1mat)
       end if

       if (any(option == [1, -2])) then
          ! Read wavefunctions.
          ! The coefficients_of_wavefunctions on file have shape [cplex, mpw, nspinor, mband, nkpt, nsppol]
         ABI_CHECK(npw == npw1, "npw != npw1 not coded")

         cg_varid = nctk_idname(wff%unwff, "coefficients_of_wavefunctions")
         ncerr = nf90_get_var(wff%unwff, cg_varid, cg(:, icg+1:), start=[1,1,1,1,ikpt,isppol], &
         count=[2,npw,nspinor,nband1,1,1])
         NCF_CHECK_MSG(ncerr, "getting cg_k")
       end if
#endif
     end if

   end if ! formeig == 1

 end if ! nband >0

!If fewer than all bands were read wind disk file forward to end of bands for this k point.
!Will have to fill the non-filled bands outside of this routine ...
 if (nband<nband_disk .or. option==-1) then
   nrec=(formeig+1)*(nband_disk-nband)
   if(option==-1)nrec=(formeig+1)*nband_disk
   call WffReadSkipRec(ios,nrec,wff)
 end if

!---------------------------------------------------------------------------
! Free memory
!---------------------------------------------------------------------------
 if (allocated(ind_cg_mpi_to_seq)) then
   ABI_DEALLOCATE(ind_cg_mpi_to_seq)
 end if

end subroutine readwf
!!***

! -------------------------------------------------------------------------------------------------

!!****f* m_rwwf/writewf
!! NAME
!! writewf
!!
!! FUNCTION
!!  This subroutine writes the block of records
!!  related to one k point, and one spin-polarization, that
!!  contains the wavefunctions (as well as the eigenvalues and occupations).
!!  The disk file unitwf should have been prepared
!!  outside of this routine, in order to write the correct records.
!!
!! INPUTS
!!  cg(2,npw*nspinor*mband)=planewave coefficients of wavefunctions,
!!  eigen((2*mband)**formeig *mband)=array for holding eigenvalues (hartree)
!!  formeig=format of the eigenvalues
!!   0 => vector of eigenvalues
!!   1 => hermitian matrix of eigenvalues
!!  icg=shift to be given to the location of the cg array
!!  ikpt=index of current k point (only needed for error message)
!!  isppol=spin polarization currently treated (only needed for error message)
!!  kg_k(3,optkg*npw)=k+g data  (only if optkg==1)
!!  mband=maximum number of bands (dimension of cg, eigen and occ)
!!  mcg=dimension of cg
!!  nband=number of bands actually in cg, eigen and occ
!!   (must be larger or equal to nband_disk, only nband_disk bands are written)
!!  nband_disk=number of bands on the disk file
!!  npw=number of plane waves
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  occ(mband)=array for holding electronic occupations
!!  option= 2 for writing cg, eigen and occ,
!!          4 for writing a file containing only eigenvalues and occupations
!!  optkg= if 1 , read or write kg_k ; if 0, do not care about kg_k
!!  wff=struct info for wavefunction
!!
!! OUTPUT
!! (none, only writing)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  WARNING : skipping k-blocks is also done in the randac subroutine
!!  WARNING : writing the two first records is also done in the dfpt_vtowfk routine
!!
!! PARENTS
!!      rwwf
!!
!! CHILDREN
!!      mpi_bcast,wffreadwrite_mpio,wffwritenpwrec,xderivewrecend
!!      xderivewrecinit,xderivewrite,xmpi_sum
!!
!! SOURCE

subroutine writewf(cg,eigen,formeig,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,&
&                  nband,nband_disk,npw,nspinor,occ,option,optkg,wff)

!Arguments ------------------------------------
 integer,intent(in) :: formeig,icg,ikpt,isppol,mband,mcg,nband,nband_disk,npw,nspinor,option,optkg
 integer,intent(in),target :: kg_k(3,optkg*npw)
 real(dp),intent(in),target :: cg(2,mcg),eigen((2*mband)**formeig*mband),occ(mband)
 type(MPI_type),intent(in) :: mpi_enreg
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
 integer :: iband,ii,ios,ipw,ispinor_index,nband2,ncid_hdr,npwso,npwsotot,npwtot,nspinortot
 integer :: unitwf,use_f90
 character(len=500) :: msg
 integer :: ikpt_this_proc,ispinor,me_cart_3d
 integer,allocatable :: ind_cg_mpi_to_seq(:)
 real(dp),ABI_CONTIGUOUS pointer :: cg_ptr(:,:)
#ifdef HAVE_NETCDF
 integer :: kg_varid,eig_varid,occ_varid,cg_varid,ncerr
 character(len=nctk_slen) :: kdep
#endif

! *********************************************************************

!Check the options
 if ( ALL(option /= (/2,4,5/)) ) then
   write(msg,'(a,i0)')' The argument option should be 2, 4 or 5. However, option=',option
   MSG_BUG(msg)
 end if

 if(wff%iomode==IO_MODE_MPI)then
   if(option==5)then
     write(msg,'(a,i0)')' With MPI-IO activated, the argument option should be 2 or 4. However, option=',option
     MSG_BUG(msg)
   end if
 end if

 if (wff%iomode==IO_MODE_ETSF) then
   ! outkss still calls this routine!
   MSG_WARNING("You should use outwff or ncwrite_cg to write data in netcdf format!")
 end if

!Check that nband_disk is not larger than nband (only for writing)
 if (nband<nband_disk) then
   write(msg,'(3a,i5,a,i5,a)') &
&   'Writing option of rwwf. One should have nband<=nband_disk',ch10,&
&   'However, nband= ',nband,', and nband_disk= ',nband_disk,'.'
   MSG_BUG(msg)
 end if

 npwtot=npw; npwso=npw*nspinor
 unitwf=wff%unwff; ncid_hdr=unitwf
 npwsotot=npwso
 nspinortot=min(2,(1+mpi_enreg%paral_spinor)*nspinor)

 use_f90=0
 if (wff%iomode==IO_MODE_FORTRAN.or.(wff%iomode ==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me)) use_f90=1

 if (wff%iomode==IO_MODE_MPI) then

   call xmpi_sum(npwsotot,wff%spaceComm_mpiio,ios)
   npwtot=npwsotot/nspinortot

   if (option/=4) then
     ABI_ALLOCATE(ind_cg_mpi_to_seq,(npwso))
     if (allocated(mpi_enreg%my_kgtab)) then
       ikpt_this_proc=mpi_enreg%my_kpttab(ikpt)
       do ispinor=1,nspinor
         ispinor_index=ispinor
         if (mpi_enreg%nproc_spinor > 1) ispinor_index = mpi_enreg%me_spinor + 1
         ind_cg_mpi_to_seq(1+npw*(ispinor-1):npw*ispinor)=npwtot*(ispinor_index-1) &
&         + mpi_enreg%my_kgtab(1:npw,ikpt_this_proc)
       end do
     else
       ind_cg_mpi_to_seq(1:npwso) = (/(ipw,ipw=1,npwso)/)
     end if
     !write(std_out,*)"MPI-IO ind_cg_mpi_to_seq", ind_cg_mpi_to_seq(1:5)
   end if
 end if

!---------------------------------------------------------------------------
!Write the first record: npw, nspinor, nband_disk
!---------------------------------------------------------------------------
!Not modified for netCDF: no need to add writing of nband_disk,npw,nspinor

 call WffWriteNpwRec(ios,nband_disk,npwtot,nspinortot,wff,opt_paral=2)

!---------------------------------------------------------------------------
!Write the second record: (k+G) vectors
!---------------------------------------------------------------------------

 if (optkg/=0.and.option/=4) then
   if(use_f90==1)then
     write(unitwf) kg_k(1:3,1:optkg*npw)
   else if (wff%iomode==IO_MODE_MPI) then

     if (allocated(mpi_enreg%my_kgtab)) then
       me_cart_3d=xmpi_comm_rank(mpi_enreg%comm_bandspinorfft)
       ikpt_this_proc=mpi_enreg%my_kpttab(ikpt)
       call xderiveWRecInit(wff,ios,me_cart_3d)
       if (mpi_enreg%me_spinor==0) then
         call xderiveWrite(wff,kg_k,3,npw,mpi_enreg%comm_bandfft, &
&         mpi_enreg%my_kgtab(1:npw,ikpt_this_proc),ios)
       end if
       call xderiveWRecEnd(wff,ios,me_cart_3d)
     else
!      MG does it work if we are not using FFT distribution ?
       call xderiveWRecInit(wff,ios )
       if (mpi_enreg%me_spinor==0) then
         call xderiveWrite(wff,kg_k,3,optkg*npw,Wff%spaceComm_mpiio,ios)
       end if
       call xderiveWRecEnd(wff,ios)
     end if

#ifdef HAVE_NETCDF
   else if (wff%iomode == IO_MODE_ETSF) then
     ! Write the reduced_coordinates_of_plane_waves for this k point.
     NCF_CHECK(nf90_inq_varid(wff%unwff, "reduced_coordinates_of_plane_waves", kg_varid))
     NCF_CHECK(nf90_get_att(wff%unwff, kg_varid, "k_dependent", kdep))
     if (kdep == "no") then
       ncerr = nf90_put_var(wff%unwff, kg_varid, kg_k, start=[1,1], count=[3,npw])
     else
       ncerr = nf90_put_var(wff%unwff, kg_varid, kg_k, start=[1,1,ikpt], count=[3,npw,1])
     end if
     NCF_CHECK_MSG(ncerr, "putting kg_k")
#endif
   end if ! end if wff%iomode
 else ! Still skip the record
   if (use_f90==1) then
     write(unitwf)
   else if (wff%iomode==IO_MODE_MPI) then
     call xderiveWRecInit(wff,wff%spaceComm_mpiio,ios)
     call xderiveWRecEnd(wff,wff%spaceComm_mpiio,ios)
   end if
 end if

!---------------------------------------------------------------------------
!Write the third record: eigenvalues and occupations
!---------------------------------------------------------------------------

!===== Case formeig=0: write eigenvalues and occupations =====
 if (formeig==0) then
   if (use_f90==1) then
     write(unitwf) (eigen(iband),iband=1,nband_disk),(occ(iband),iband=1,nband_disk)
   else if(wff%iomode==IO_MODE_MPI) then
     if (wff%me_mpiio==0) then
       call xderiveWRecInit(wff,ios)
       call xderiveWrite(wff,eigen,nband_disk,xmpi_comm_self,ios)
       call xderiveWrite(wff,occ,nband_disk,xmpi_comm_self,ios)
       call xderiveWRecEnd(wff,ios)
     end if
#ifdef HAVE_MPI
     call MPI_BCAST(wff%offwff,1,wff%offset_mpi_type,0,wff%spaceComm_mpiio,ios)
#endif
#ifdef HAVE_NETCDF
   else if (wff%iomode == IO_MODE_ETSF) then
     ! Write eigenvalues and occupation factors.
     NCF_CHECK(nf90_inq_varid(wff%unwff, "eigenvalues", eig_varid))
     ncerr = nf90_put_var(wff%unwff, eig_varid, eigen, start=[1,ikpt,isppol], count=[mband,1,1])
     NCF_CHECK_MSG(ncerr, "putting eig_k")

     NCF_CHECK(nf90_inq_varid(wff%unwff, "occupations", occ_varid))
     ncerr = nf90_put_var(wff%unwff, occ_varid, occ, start=[1,ikpt,isppol], count=[mband,1,1])
     NCF_CHECK_MSG(ncerr, "putting occ_k")
#endif
   end if

!  ===== Case formeig=1: write matrix of eigenvalues =====
!  Will be written later (together with wave-functions)
 else if(formeig==1)then
 end if ! formeig

!---------------------------------------------------------------------------
!Write the wave-function coefficients
!---------------------------------------------------------------------------

!===== Case formeig=0: write only wave-functions =====
 if (formeig==0) then
!  If option=4, do not write wave functions
   if (option/=4) then
     if (use_f90==1) then
       do iband=1,nband_disk
         ipw=(iband-1)*npwso+icg
         if(option/=5)then
           write(unitwf) cg(1:2,ipw+1:ipw+npwso) ! VALGRIND complains some elements of cg are not initialized, but written
         else
           write(unitwf)
         end if
       end do
     else if(wff%iomode==IO_MODE_MPI)then
       cg_ptr => cg ! Need pointer to bypass "inout" intent attribute
       call WffReadWrite_mpio(wff,2,cg_ptr,mcg,icg,nband_disk,npwso,npwsotot,ind_cg_mpi_to_seq,ios)
       nullify(cg_ptr)
#ifdef HAVE_NETCDF
     else if (wff%iomode == IO_MODE_ETSF .and. option/=5) then

       NCF_CHECK(nf90_inq_varid(wff%unwff, "reduced_coordinates_of_plane_waves", kg_varid))
       NCF_CHECK(nf90_get_att(wff%unwff, kg_varid, "k_dependent", kdep))

       ! The coefficients_of_wavefunctions on file have shape [cplex, mpw, nspinor, mband, nkpt, nsppol]
       NCF_CHECK(nf90_inq_varid(wff%unwff, "coefficients_of_wavefunctions", cg_varid))
       !write(std_out,*)"put cg, count: ",[2,npw,nspinor,nband,1,1]
       ncerr = nf90_put_var(wff%unwff, cg_varid, cg(:, icg+1:), start=[1,1,1,1,ikpt,isppol], &
       count=[2,npw,nspinor,nband,1,1])
       NCF_CHECK_MSG(ncerr, "putting cg_k")
       !write(std_out,*)"after cg"
#endif
     end if
   end if ! option/=4

!  ===== Case formeig=1: write eigenvalues and wave-functions =====
 else if(formeig==1)then

!  Not available for NETCDF and ETSF_IO
   ABI_CHECK(wff%iomode /= IO_MODE_ETSF, "ETSF-write-eigen1 not coded!")

!  ABI_CHECK(nband_disk==nband,"nband_disk!=nband")

   nband2=2*nband_disk
   do iband=1,nband_disk
     ipw=(iband-1)*npwso+icg
     ii=(iband-1)*nband2
     if(use_f90==1)then
       write(unitwf) eigen(1+ii:nband2+ii)
       if (option/=5) then
         write(unitwf) cg(1:2,1+ipw:npwso+ipw)
       else
         write(unitwf)
       end if
     else if(wff%iomode==IO_MODE_MPI)then
!      Should use an access with a "view"
       call xderiveWRecInit(wff,ios)
       call xderiveWrite(wff,eigen(1+ii:ii+nband2),nband2,wff%spaceComm_mpiio,ios)
       call xderiveWRecEnd(wff,ios)
       if (option/=4) then
         call xderiveWRecInit(wff,ios)
         call xderiveWrite(wff,cg(1:2,ipw+1:ipw+npwso),2,npwso,wff%spaceComm_mpiio,ios)
         call xderiveWRecEnd(wff,ios)
       end if
     end if
   end do

 end if ! formeig

!---------------------------------------------------------------------------
!Final statements
!---------------------------------------------------------------------------

 if (allocated(ind_cg_mpi_to_seq)) then
   ABI_DEALLOCATE(ind_cg_mpi_to_seq)
 end if

 RETURN

!Silence compiler warning.
 ABI_UNUSED((/ii,mpi_enreg%me/))

end subroutine writewf
!!***

end module m_rwwf
!!***
