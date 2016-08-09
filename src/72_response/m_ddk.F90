!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ddk
!! NAME
!!  m_ddk
!!
!! FUNCTION
!!  Objects and methods to extract data from DDK files.
!!  The DDK files are binary (soon also netcdf) files with Hamiltonian derivatives
!!  wrt k, and the corresponding wave functions
!!
!! COPYRIGHT
!! Copyright (C) 2016 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_ddk

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_distribfft
 use m_nctk
#ifdef HAVE_TRIO_NETCDF
 use netcdf
#endif
 use m_hdr
 use m_kptrank
 use m_fstab
 use m_wfk

 use m_fstrings,      only : sjoin, itoa
 use m_io_tools,      only : open_file, file_exists, iomode_from_fname
! use m_numeric_tools, only : wrap2_pmhalf, vdiff_eval, vdiff_print
! use m_copy,          only : alloc_copy
 use defs_abitypes,   only : hdr_type, mpi_type
 use m_mpinfo,        only : destroy_mpi_enreg
 use m_crystal,       only : crystal_t, crystal_free
 use m_crystal_io,    only : crystal_from_hdr

 implicit none

 private
!!***

 integer,public,parameter :: ddk_last_version = 1

 integer,private,parameter :: DDK_NOMODE    = 0
 integer,private,parameter :: DDK_READMODE  = 1
 integer,private,parameter :: DDK_WRITEMODE = 2

 ! FIXME
 !real(dp),public,parameter :: DDB_QTOL=2.0d-8
 ! Tolerance for the identification of two wavevectors

 integer,private,parameter :: FPOS_EOF = -1



!----------------------------------------------------------------------

!!****t* m_ddk/ddk_t
!! NAME
!!  ddk_t
!!
!! FUNCTION
!!  object containing ddk derivatives ([H,r] proportional to band velocities)
!!
!! NOTES
!!
!! SOURCE

 type,public :: ddk_t

  integer :: fh(3)
  ! file handler
  !  Fortran unit number if iomode==IO_MODE_FORTRAN
  !  MPI file handler if iomode==IO_MODE_MPI
  ! TODO: add netcdf?

  integer :: comm
  ! MPI communicator used for IO.

  integer :: version
  ! File format version read from file.

  integer :: iomode=IO_MODE_FORTRAN
  ! Method used to access the DDK file:
  !   IO_MODE_FORTRAN for usual Fortran IO routines
  !   IO_MODE_MPI if MPI/IO routines.
  !   IO_MODE_ETSF netcdf files in etsf format

  integer :: rw_mode = DDK_NOMODE
   ! (Read|Write) mode 

  integer :: current_fpos
  ! The current position of the file pointer used for sequential access with Fortran-IO
  !  FPOS_EOF signals the end of file

  integer :: nsppol
   ! Number of spin polarizations.

  integer :: nspinor
   ! Number of spinor components.

  integer :: nkfs
    ! Number of k-points on the Fermi-surface (FS-BZ).

  integer :: maxnb
   ! Max number of bands on the FS.

  integer :: usepaw
   ! 1 if PAW calculation, 0 otherwise

  integer :: prtvol=0
   ! Verbosity level

  logical :: debug=.False.
   ! Debug flag

  character(len=fnlen) :: path(3) = ABI_NOFILE
   ! File name 

  real(dp) :: acell(3),rprim(3,3),gprim(3,3)

  real(dp), allocatable :: velocity (:,:,:,:)
   ! dims (3,maxnb,nkfs,nsppol)
 
   ! Used for slow FT.
  type(crystal_t) :: cryst

  type(mpi_type) :: mpi_enreg

 end type ddk_t

 public :: ddk_init              ! Initialize the object.
 public :: ddk_read_from_file    ! Read ddk matrix elements from file
 public :: ddk_free              ! Close the file and release the memory allocated.
 public :: ddk_print             ! output values

!!***

CONTAINS
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddk_init
!! NAME
!!  ddk_init
!!
!! FUNCTION
!!  Initialize the object from file. This is a COLLECTIVE procedure that must be called 
!!  by each process in the MPI communicator comm.
!!
!! INPUTS
!!   path=Filename.
!!   comm=MPI communicator.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddk_init(ddk, path, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddk_init'
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: path(3)
 integer,intent(in) :: comm
 type(ddk_t),intent(inout) :: ddk

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,timrev2=2
 integer :: ierr,fform
 integer :: fform_ref
 integer :: my_rank
 integer :: restart, restartpaw
 character(len=500) :: msg
 type(hdr_type) :: hdr1,hdr_ref
!arrays 
 real(dp), allocatable :: eigen1(:)

!************************************************************************

 my_rank = xmpi_comm_rank(comm) 
 ddk%path = path; ddk%comm = comm
 ddk%iomode = iomode_from_fname(path(1))

! if (ddk%iomode == IO_MODE_MPI) then
!   write (msg,'(a)') " MPI iomode found for ddk - not coded yet, will default to FORTRAN "
!   MSG_WARNING(msg)
!   ddk%iomode = IO_MODE_FORTRAN
! end if


! in this call everything is broadcast properly to the whole comm
 call wfk_read_h1mat(path(1), eigen1, hdr_ref, comm)
 fform_ref = 2
 ABI_DEALLOCATE(eigen1)

 if (ddk%debug) call hdr_echo(hdr_ref,fform_ref,4,unit=std_out)

 ! check that the other 2 headers are compatible
 call wfk_read_h1mat(path(2), eigen1, hdr1, comm)
 fform = 2
 ABI_DEALLOCATE(eigen1)
 call hdr_check(fform,fform_ref,hdr1,hdr_ref,'COLL',restart,restartpaw)
 call hdr_free(hdr1)


 call wfk_read_h1mat(path(3), eigen1, hdr1, comm)
 fform = 2
 ABI_DEALLOCATE(eigen1)
 call hdr_check(fform,fform_ref,hdr1,hdr_ref,'COLL',restart,restartpaw)
 call hdr_free(hdr1)


 ddk%nsppol = hdr_ref%nsppol
 ddk%nspinor = hdr_ref%nspinor
 ddk%usepaw = hdr_ref%usepaw
 ABI_CHECK(ddk%usepaw == 0, "PAW not yet supported")

 ! Init crystal_t
 call crystal_from_hdr(ddk%cryst,hdr_ref,timrev2)
 call hdr_free(hdr_ref)

 ! Compute rprim, and gprimd. Used for slow FFT q--r if multiple shifts
 call mkradim(ddk%acell,ddk%rprim,ddk%cryst%rprimd)
 call matr3inv(ddk%rprim,ddk%gprim)

 ! MPI_type needed for calling fourdp!
 call initmpi_seq(ddk%mpi_enreg)

end subroutine ddk_init
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddk_read_from_file
!! NAME
!!  ddk_read_from_file
!!
!! FUNCTION
!!  Read in ddk matrix elements from files
!!
!! INPUTS
!!   comm=MPI communicator
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddk_read_from_file(comm, ddk, fstab)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddk_read_from_file'
 use interfaces_14_hidewrite
 use interfaces_72_response
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 type(ddk_t),intent(inout) :: ddk
 type(fstab_t),target,intent(in) :: fstab(ddk%nsppol)

!arrays 

!Local variables-------------------------------
!scalars
 integer :: nprocs, my_rank, master=0
 integer :: idir
 integer :: ierr
 integer :: ikfs, ikpt, isppol, ik_ibz
 integer :: symrankkpt, ikpt_ddk
 integer :: iband, bd2tot_index
 integer :: mband
 integer :: bstart_k, nband_k
 integer :: nband_in
 integer :: vdim
 character(len=500) :: message
 real(dp) :: tmpveloc(3), tmpveloc2(3)
 type(hdr_type) :: hdr1
 real(dp), allocatable :: eigen1(:)
 real(dp), allocatable :: velocityp(:,:)
 real(dp), allocatable :: velocitypp(:,:)
 type(kptrank_type) :: kptrank_t
 type(fstab_t), pointer :: fs

!************************************************************************

 if (ddk%rw_mode /= ddk_NOMODE) then
   MSG_ERROR("ddk should be in ddk_NOMODE before open_read is called.")
 end if
 ddk%rw_mode = ddk_READMODE

 nprocs = xmpi_comm_size(comm)
 my_rank = xmpi_comm_rank(comm)

 ddk%maxnb = maxval(fstab(:)%maxnb)
 ddk%nkfs = maxval(fstab(:)%nkfs)
 ABI_ALLOCATE(ddk%velocity, (3,ddk%maxnb,ddk%nkfs,ddk%nsppol))

! TODO: make this parallel enabled - at least check for master, at best get distributed k

 do idir = 1,3
   ! Open the files
   select case (ddk%iomode)
   case (IO_MODE_FORTRAN)
#ifdef DEV_MG_WFK
   write(message,'(a,i2)')' NEW DDK FILES, iomode = ',ddk%iomode
   call wrtout(std_out,message,'COLL')
   call wfk_read_h1mat (ddk%path(idir), eigen1, hdr1, comm)
#else 
   write(message,'(a,i2)')' OLD DDK FILES, iomode = ',ddk%iomode
   call wrtout(std_out,message,'COLL')
   if (my_rank == master) then
     call inpgkk(eigen1,ddk%path(idir),hdr1)
   end if
   if (xmpi_comm_size(comm) > 1) then
       call hdr_bcast(hdr1, master, my_rank, comm)
       if (my_rank /= master) then
         mband = maxval(hdr1%nband)
         ABI_ALLOCATE(eigen1, (2*hdr1%nsppol*hdr1%nkpt*mband**2))
       end if
       call xmpi_bcast(eigen1, master, comm, ierr)
   end if
#endif 
! TODO: later combine these case statements, which can call the same routine, and eliminate the old ddk file format
   case (IO_MODE_ETSF, IO_MODE_MPI)
     write(message,'(a,i2)')' NEW DDK FILES, iomode = ',ddk%iomode
     call wrtout(std_out,message,'COLL')
     call wfk_read_h1mat (ddk%path(idir), eigen1, hdr1, comm)
!   case (IO_MODE_MPI)
!     MSG_ERROR("MPI not coded")
   case default
     MSG_ERROR(sjoin("Unsupported iomode:", itoa(ddk%iomode)))
   end select
  
   nband_in = maxval(hdr1%nband)

!need correspondence hash between the DDK and the fs k-points
   call mkkptrank (hdr1%kptns,hdr1%nkpt,kptrank_t)
   do isppol=1,ddk%nsppol
     fs => fstab(isppol)

     do ikfs=1,fs%nkfs
       ik_ibz = fs%istg0(1,ikfs)
       call get_rank_1kpt (fs%kpts(:,ikfs),symrankkpt, kptrank_t)
       ikpt_ddk = kptrank_t%invrank(symrankkpt)
       if (ikpt_ddk == -1) then
         write(std_out,*)'ddk_read_from_file ******** error in correspondence between ddk and fsk kpoint sets'
         write(std_out,*)' kpt sets in fsk and ddk files must agree.'
         MSG_ERROR("Aborting now")
       end if

       bd2tot_index=2*nband_in**2*(ikpt_ddk-1) + 2*nband_in**2*hdr1%nkpt*(isppol-1)
       bstart_k = fs%bstcnt_ibz(1, ik_ibz)
       nband_k = fs%bstcnt_ibz(2, ik_ibz)
!      first derivative eigenvalues for k-point. Diagonal of eigen1 is real -> only use that part
       do iband = bstart_k, bstart_k+nband_k-1
         ddk%velocity(idir, iband-bstart_k+1, ikfs, isppol)=eigen1(bd2tot_index + 2*nband_in*(iband-1) + iband)
       end do

     end do
   end do
   call destroy_kptrank (kptrank_t)


! This should always be allocated (from within the read h1mat or inpgkk subroutines)
!   if (allocated(eigen1)) then
     ABI_DEALLOCATE(eigen1)
!   end if
   call hdr_free(hdr1)
 end do

 ! process the eigenvalues(1): rotate to cartesian and divide by 2 pi
 ! use DGEMM better here on whole matrix and then reshape?
 vdim = ddk%maxnb*ddk%nkfs*ddk%nsppol
 ABI_ALLOCATE(velocityp, (3,vdim))
 velocityp = reshape (ddk%velocity, (/3,vdim/))
 ABI_ALLOCATE(velocitypp, (3,vdim))
 velocitypp = zero
 call dgemm('n','n',3,vdim,3,one,ddk%cryst%rprimd,3,velocityp,3,zero,velocitypp,3)
 ABI_DEALLOCATE(velocityp)

! do isppol = 1, ddk%nsppol
!   do ikpt = 1, ddk%nkfs
!     do iband = 1, ddk%maxnb
!       tmpveloc = ddk%velocity (:, iband, ikpt, isppol)
!       call dgemm('n','n',3,3,3,one,ddk%cryst%rprimd,3,tmpveloc,3,zero,tmpveloc2,3)
!       ddk%velocity (:, iband, ikpt, isppol) = tmpveloc2
!     end do
!   end do
! end do

 ddk%velocity = reshape (velocitypp, (/3,ddk%maxnb,ddk%nkfs,ddk%nsppol/))
 ABI_DEALLOCATE(velocitypp)

 call destroy_kptrank (kptrank_t)

end subroutine ddk_read_from_file
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddk_free
!! NAME
!!  ddk_free
!!
!! FUNCTION
!! Close the file and release the memory allocated.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddk_free(ddk)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddk_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ddk_t),intent(inout) :: ddk

!************************************************************************

 ! integer arrays

 ! real arrays
 if (allocated(ddk%velocity)) then
   ABI_FREE(ddk%velocity)
 end if
 
 ! types
 call crystal_free(ddk%cryst)
 call destroy_mpi_enreg(ddk%mpi_enreg)

end subroutine ddk_free
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddk_print
!! NAME
!!  ddk_print
!!
!! FUNCTION
!!  Print info on the object.
!!
!! INPUTS
!! [unit]=the unit number for output 
!! [prtvol]=verbosity level
!! [mode_paral]=either "COLL" or "PERS"
!!
!! OUTPUT
!!  Only printing.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddk_print(ddk, header, unit, prtvol, mode_paral)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddk_print'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol,unit
 character(len=4),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header
 type(ddk_t),intent(in) :: ddk

!Local variables-------------------------------
!scalars
 integer :: my_unt,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_unt =std_out; if (PRESENT(unit)) my_unt   =unit
 my_prtvol=0    ; if (PRESENT(prtvol)) my_prtvol=prtvol
 my_mode='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral
                                                                    
 msg=' ==== Info on the ddk% object ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 write(std_out,*)"Number of FS bands: ",ddk%maxnb
 write(std_out,*)"Number of FS k-points: ",ddk%nkfs
 write(std_out,*)"Number of spin channels: ",ddk%nsppol
 write(std_out,*)"Path to files: "
 write(std_out,*)"  ", trim(ddk%path(1))
 write(std_out,*)"  ", trim(ddk%path(2))
 write(std_out,*)"  ", trim(ddk%path(3))

end subroutine ddk_print
!!***


END MODULE m_ddk
