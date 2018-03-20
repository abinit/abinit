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
!! Copyright (C) 2016-2018 ABINIT group (MJV)
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
 use m_nctk
 use m_hdr
 use m_kptrank
 use m_fstab
 use m_wfk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_fstrings,      only : sjoin, itoa, endswith
 use m_io_tools,      only : iomode_from_fname
 use defs_abitypes,   only : hdr_type
 use m_crystal,       only : crystal_t, crystal_free
 use m_crystal_io,    only : crystal_from_hdr
 use defs_datatypes,  only : ebands_t

 implicit none

 private
!!***

 integer,private,parameter :: DDK_NOMODE    = 0
 integer,private,parameter :: DDK_READMODE  = 1
 integer,private,parameter :: DDK_WRITEMODE = 2

!----------------------------------------------------------------------

!!****t* m_ddk/ddk_t
!! NAME
!!  ddk_t
!!
!! FUNCTION
!!  object containing ddk derivatives ([H,r] proportional to band velocities)
!!
!! SOURCE

 type,public :: ddk_t

  !integer :: fh(3)
  ! file handler
  !  Fortran unit number if iomode==IO_MODE_FORTRAN
  !  MPI file handler if iomode==IO_MODE_MPI

  integer :: comm
  ! MPI communicator used for IO.

  !integer :: version
  ! File format version read from file.

  integer :: iomode=IO_MODE_FORTRAN
  ! Method used to access the DDK file:
  !   IO_MODE_FORTRAN for usual Fortran IO routines
  !   IO_MODE_MPI if MPI/IO routines.
  !   IO_MODE_ETSF netcdf files in etsf format

  integer :: rw_mode = DDK_NOMODE
   ! (Read|Write) mode

  !integer :: current_fpos
  ! The current position of the file pointer used for sequential access with Fortran-IO
  !  FPOS_EOF signals the end of file

  integer :: nsppol
   ! Number of spin polarizations.

  integer :: nspinor
   ! Number of spinor components.

  integer :: nene
    ! Number of energy points we may need the fsavg at

  integer :: nkfs
    ! Number of k-points on the Fermi-surface (FS-BZ).

  integer :: maxnb
   ! Max number of bands on the FS.
   ! TODO: Maybe maxnbfs

  integer :: usepaw
   ! 1 if PAW calculation, 0 otherwise

  integer :: prtvol=0
   ! Verbosity level

  logical :: debug=.False.
   ! Debug flag

  character(len=fnlen) :: paths(3) = ABI_NOFILE
   ! File name

  real(dp) :: acell(3),rprim(3,3),gprim(3,3)
   ! TODO: Are these really needed?

  real(dp), allocatable :: velocity (:,:,:,:)
   ! (3,maxnb,nkfs,nsppol)
   ! velocity on the FS in cartesian coordinates.

  real(dp), allocatable :: velocity_fsavg (:,:,:)
   ! (3,nene,nsppol)
   ! velocity on the FS in cartesian coordinates.

  logical :: use_ncddk(3)
   ! True if we are readin DDK matrix elements from DDK.nc instead of WFK file

  type(crystal_t) :: cryst
   ! Crystal structure read from file

 end type ddk_t

 public :: ddk_init              ! Initialize the object.
 public :: ddk_read_fsvelocities ! Read FS velocities from file.
 public :: ddk_fs_average_veloc  ! find FS average of velocity squared
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
!!   paths=3 Filenames (could be either DFPT WFK files of DKK.nc files.
!!   comm=MPI communicator.
!!
!! PARENTS
!!      eph
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine ddk_init(ddk, paths, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddk_init'
 use interfaces_32_util
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: paths(3)
 integer,intent(in) :: comm
 type(ddk_t),intent(inout) :: ddk

!Local variables-------------------------------
!scalars
 integer,parameter :: timrev2=2,fform2=2
 integer :: my_rank, restart, restartpaw, ii
!arrays
 integer :: fforms(3)
 type(hdr_type) :: hdrs(3)

!************************************************************************

 my_rank = xmpi_comm_rank(comm)
 ddk%paths = paths; ddk%comm = comm
 ddk%iomode = iomode_from_fname(paths(1))

 ! In this calls everything is broadcast properly to the whole comm
 do ii=1,3
   ddk%use_ncddk(ii) = endswith(paths(ii), "_DDK.nc")
   call hdr_read_from_fname(hdrs(ii), paths(ii), fforms(ii), comm)
   if (ddk%debug) call hdr_echo(hdrs(ii), fforms(ii), 4, unit=std_out)
   ! check that 2 headers are compatible
   if (ii > 1) call hdr_check(fform2, fform2, hdrs(ii-1), hdrs(ii), 'COLL', restart, restartpaw)
 end do

 ddk%nsppol = hdrs(1)%nsppol
 ddk%nspinor = hdrs(1)%nspinor
 ddk%usepaw = hdrs(1)%usepaw
 ABI_CHECK(ddk%usepaw == 0, "PAW not yet supported")

 ! Init crystal_t
 call crystal_from_hdr(ddk%cryst, hdrs(1), timrev2)

 ! Compute rprim, and gprimd. Used for slow FFT q--r if multiple shifts
 call mkradim(ddk%acell,ddk%rprim,ddk%cryst%rprimd)
 call matr3inv(ddk%rprim,ddk%gprim)

 do ii=1,3
   call hdr_free(hdrs(ii))
 end do

end subroutine ddk_init
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddk_read_fsvelocities
!! NAME
!!  ddk_read_fsvelocities
!!
!! FUNCTION
!!  Read FS velocities from DDK files. Returned in ddk%velocity
!!
!! INPUTS
!!   fstab(ddk%nsppol)=Tables with the correspondence between points of the Fermi surface (FS)
!!     and the k-points in the IBZ
!!   comm=MPI communicator
!!
!! PARENTS
!!      m_phgamma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine ddk_read_fsvelocities(ddk, fstab, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddk_read_fsvelocities'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 type(ddk_t),intent(inout) :: ddk
 type(fstab_t),target,intent(in) :: fstab(ddk%nsppol)

!Local variables-------------------------------
!scalars
 integer :: idir, ikfs, isppol, ik_ibz, ii
 integer :: symrankkpt, ikpt_ddk,iband, bd2tot_index
 integer :: bstart_k, nband_k, nband_in, vdim
#ifdef HAVE_NETCDF
 integer :: ncid, varid, nc_fform
#endif
 type(hdr_type) :: hdr1
 type(kptrank_type) :: kptrank_t
 type(fstab_t), pointer :: fs
 character(len=500) :: msg
!arrays
 real(dp), allocatable :: eigen1(:)
 real(dp), allocatable :: velocityp(:,:)

!************************************************************************

 if (ddk%rw_mode /= ddk_NOMODE) then
   MSG_ERROR("ddk should be in ddk_NOMODE before open_read is called.")
 end if
 ddk%rw_mode = ddk_READMODE

 ddk%maxnb = maxval(fstab(:)%maxnb)
 ddk%nkfs = maxval(fstab(:)%nkfs)
 ABI_MALLOC(ddk%velocity, (3,ddk%maxnb,ddk%nkfs,ddk%nsppol))

 call wrtout(std_out, sjoin('Read DDK FILES with iomode=', itoa(ddk%iomode)), 'COLL')
 do idir = 1,3
   ! Open the files. All procs in comm receive hdr1 and eigen1
   if (ddk%use_ncddk(idir)) then
#ifdef HAVE_NETCDF
     NCF_CHECK(nctk_open_read(ncid, ddk%paths(ii), comm))
     call hdr_ncread(hdr1, ncid, nc_fform)
     varid = nctk_idname(ncid, "h1_matrix_elements")
     ABI_MALLOC(eigen1, (2*hdr1%mband*hdr1%mband*hdr1%nkpt*ddk%nsppol))
     !NCF_CHECK(nctk_set_collective(ncid, varid))
     NCF_CHECK(nf90_get_var(ncid, varid, eigen1, count=[2, hdr1%mband, hdr1%mband, hdr1%nkpt, hdr1%nsppol]))
     NCF_CHECK(nf90_close(ncid))
#else
     MSG_ERROR("Netcdf not available!")
#endif
   else
     call wfk_read_h1mat(ddk%paths(idir), eigen1, hdr1, comm)
   end if
   nband_in = maxval(hdr1%nband)

   ! need correspondence hash between the DDK and the fs k-points
   call mkkptrank (hdr1%kptns,hdr1%nkpt,kptrank_t)
   do isppol=1,ddk%nsppol
     fs => fstab(isppol)
     do ikfs=1,fs%nkfs
       ik_ibz = fs%istg0(1,ikfs)
       call get_rank_1kpt (fs%kpts(:,ikfs),symrankkpt, kptrank_t)
       ikpt_ddk = kptrank_t%invrank(symrankkpt)
       if (ikpt_ddk == -1) then
         write(msg, "(3a)")&
           "Error in correspondence between ddk and fsk kpoint sets",ch10,&
           "kpt sets in fsk and ddk files must agree."
         MSG_ERROR(msg)
       end if

       bd2tot_index=2*nband_in**2*(ikpt_ddk-1) + 2*nband_in**2*hdr1%nkpt*(isppol-1)
       bstart_k = fs%bstcnt_ibz(1, ik_ibz)
       nband_k = fs%bstcnt_ibz(2, ik_ibz)
       ! first derivative eigenvalues for k-point. Diagonal of eigen1 is real -> only use that part
       do iband = bstart_k, bstart_k+nband_k-1
         ddk%velocity(idir, iband-bstart_k+1, ikfs, isppol)=eigen1(bd2tot_index + 2*nband_in*(iband-1) + iband)
       end do
     end do
   end do

   ABI_FREE(eigen1)
   call destroy_kptrank(kptrank_t)
   call hdr_free(hdr1)
 end do ! idir

 ! process the eigenvalues(1): rotate to cartesian and divide by 2 pi
 ! use DGEMM better here on whole matrix and then reshape?
 vdim = ddk%maxnb*ddk%nkfs*ddk%nsppol
 ABI_MALLOC(velocityp, (3,vdim))
 velocityp = zero
 call dgemm('n','n',3,vdim,3,one,ddk%cryst%rprimd,3,ddk%velocity,3,zero,velocityp,3)

! do isppol = 1, ddk%nsppol
!   do ikpt = 1, ddk%nkfs
!     do iband = 1, ddk%maxnb
!       tmpveloc = ddk%velocity (:, iband, ikpt, isppol)
!       call dgemm('n','n',3,3,3,one,ddk%cryst%rprimd,3,tmpveloc,3,zero,tmpveloc2,3)
!       ddk%velocity (:, iband, ikpt, isppol) = tmpveloc2
!     end do
!   end do
! end do

 ddk%velocity = reshape (velocityp, [3,ddk%maxnb,ddk%nkfs,ddk%nsppol])

 ABI_FREE(velocityp)
 call destroy_kptrank (kptrank_t)

end subroutine ddk_read_fsvelocities
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddk_fs_average_veloc
!! NAME
!!  ddk_fs_average_veloc
!!
!! FUNCTION
!!  Perform Fermi surface average of velocity squared then square rooted, print and store in ddk object
!!
!! INPUTS
!!   ddk = object with electron band velocities
!!   fstab(ddk%nsppol)=Tables with the correspondence between points of the Fermi surface (FS)
!!     and the k-points in the IBZ
!!   comm=MPI communicator
!!
!! PARENTS
!!      m_phgamma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine ddk_fs_average_veloc(ddk, ebands, fstab, sigmas)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddk_fs_average_veloc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!integer,intent(in) :: comm  ! could distribute this over k in the future
 real(dp),intent(in) :: sigmas(:)
 type(ebands_t),intent(in) :: ebands
 type(ddk_t),intent(inout) :: ddk
 type(fstab_t),target,intent(in) :: fstab(ddk%nsppol)

!Local variables-------------------------------
!scalars
 integer :: idir, ikfs, isppol, ik_ibz, iene
 integer :: iband
 integer :: mnb, nband_k
 integer :: nsig
 type(fstab_t), pointer :: fs
!arrays
 real(dp), allocatable :: wtk(:,:)

!************************************************************************

 ddk%nene = fstab(1)%nene
 ABI_MALLOC(ddk%velocity_fsavg, (3,ddk%nene,ddk%nsppol))
 ddk%velocity_fsavg = zero

 nsig = size(sigmas, dim=1)
 mnb = 1
 do isppol=1,ddk%nsppol
   fs => fstab(isppol)
   mnb = max(mnb, maxval(fs%bstcnt_ibz(2, :)))
 end do
 ABI_MALLOC(wtk, (nsig,mnb))

 do isppol=1,ddk%nsppol
   fs => fstab(isppol)
   do iene = 1, fs%nene
     do ikfs=1,fs%nkfs
       ik_ibz = fs%istg0(1,ikfs)
       nband_k = fs%bstcnt_ibz(2, ik_ibz)
       call fstab_weights_ibz(fs, ebands, ik_ibz, isppol, sigmas, wtk, iene)

       do idir = 1,3
         do iband = 1, nband_k
           ddk%velocity_fsavg(idir, iene, isppol) = ddk%velocity_fsavg(idir, iene, isppol) + &
&             wtk(1,iband) * ddk%velocity(idir, iband, ikfs, isppol)**2
!&             fs%tetra_wtk_ene(iband,ik_ibz,iene) * ddk%velocity(idir, iband, ikfs, isppol)**2
         end do
       end do ! idir
     end do ! ikfs
   end do ! iene
  ! sqrt is element wise on purpose
   ddk%velocity_fsavg(:,:,isppol) = sqrt(ddk%velocity_fsavg(:,:,isppol)) / dble(fs%nkfs)
 end do ! isppol

 ABI_DEALLOCATE(wtk)

end subroutine ddk_fs_average_veloc
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
!!      eph
!!
!! CHILDREN
!!      wrtout
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
   ABI_DEALLOCATE(ddk%velocity)
 end if
 if (allocated(ddk%velocity_fsavg)) then
   ABI_DEALLOCATE(ddk%velocity_fsavg)
 end if

 ! types
 call crystal_free(ddk%cryst)

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
!!      wrtout
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
 write(std_out,*)"Paths to files: "
 write(std_out,*)"  ", trim(ddk%paths(1))
 write(std_out,*)"  ", trim(ddk%paths(2))
 write(std_out,*)"  ", trim(ddk%paths(3))

end subroutine ddk_print
!!***

END MODULE m_ddk
