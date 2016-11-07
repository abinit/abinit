!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dvdb
!! NAME
!!  m_dvdb
!!
!! FUNCTION
!!  Objects and methods to extract data from the DVDB file.
!!  The DVDB file is binary file with a collection of DFPT potentials
!!  associated to the different perturbations (idir, ipert, qpt).
!!  DVDB files are produced with the `mrgdv` utility.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2016 ABINIT group (MG)
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

MODULE m_dvdb

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_distribfft
 use m_nctk
 use m_sort
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_hdr

 use m_fstrings,      only : strcat, sjoin, itoa, ktoa, endswith
 use m_io_tools,      only : open_file, file_exists
 use m_numeric_tools, only : wrap2_pmhalf, vdiff_eval, vdiff_print
 use m_copy,          only : alloc_copy
 use defs_abitypes,   only : hdr_type, mpi_type
 use m_mpinfo,        only : destroy_mpi_enreg
 use m_fftcore,       only : ngfft_seq
 use m_fft_mesh,      only : rotate_fft_mesh, times_eigr, times_eikr, ig2gfft
 use m_crystal,       only : crystal_t, crystal_free
 use m_crystal_io,    only : crystal_from_hdr

 implicit none

 private
!!***

 integer,public,parameter :: dvdb_last_version = 1

 integer,private,parameter :: DVDB_NOMODE    = 0
 integer,private,parameter :: DVDB_READMODE  = 1
 integer,private,parameter :: DVDB_WRITEMODE = 2

 ! FIXME
 !real(dp),public,parameter :: DDB_QTOL=2.0d-8
 ! Tolerance for the identification of two wavevectors

 integer,private,parameter :: FPOS_EOF = -1

!----------------------------------------------------------------------

!!****t* m_dvdb/dvdb_t
!! NAME
!!  dvdb_t
!!
!! FUNCTION
!!  Database of DFPT results. The database contains `numv1` perturbations
!!  and the corresponding first order local potential in real space on the FFT mesh.
!!  Note that one can have different FFT meshes for the different perturbations
!!
!! NOTES
!!  natom, nspden, nspinor, and usepaw are global variables in the sense that it's not possible to add
!!  new entries to the database if these dimension differ from the global ones.
!!
!! SOURCE

 type,public :: dvdb_t

  integer :: fh
   ! file handler
   !  Fortran unit number if iomode==IO_MODE_FORTRAN
   !  MPI file handler if iomode==IO_MODE_MPI

  integer :: comm
  ! MPI communicator used for IO.

  integer :: version
  ! File format version read from file.

  integer :: iomode=IO_MODE_FORTRAN
  ! Method used to access the DVDB file:
  !   IO_MODE_FORTRAN for usual Fortran IO routines
  !   IO_MODE_MPI if MPI/IO routines.

  integer :: rw_mode = DVDB_NOMODE
   ! (Read|Write) mode

  integer :: current_fpos
  ! The current position of the file pointer used for sequential access with Fortran-IO
  !  FPOS_EOF signals the end of file

  integer :: numv1
  ! Total number of v1 potentials present in file.

  integer :: nqpt
  ! Number of q-points

  integer :: natom
   ! Number of atoms

  integer :: natom3
   ! 3 * natom

  integer :: nspden
   ! Number of spin density components

  integer :: nsppol
   ! Number of spin polarizations.

  integer :: nspinor
   ! Number of spinor components.

  integer :: usepaw
   ! 1 if PAW calculation, 0 otherwise

  integer :: mpert
   ! Maximum number of perturbations

  integer :: prtvol=0
   ! Verbosity level

  logical :: debug=.False.
   ! Debug flag

  logical :: symv1=.False.
   ! Activate symmetrization of v1 potentials.

  character(len=fnlen) :: path = ABI_NOFILE
   ! File name

  integer,allocatable :: pos_dpq(:,:,:)
   ! pos_dpq(3, mpert, nqpt)
   ! The position of the (idir, ipert, iqpt) potential in the file
   ! 0 if the corresping entry is not available.

  integer,allocatable :: cplex_v1(:)
  ! cplex_v1(numv1)
  ! The value of cplex for each v1(cplex*nfft, nspden) potential
  ! 2 if the potential is complex, 1 if real (q==Gamma)

  integer :: ngfft(18)=-1
   ! Info on the FFT to be used for the potentials.

  ! FIXME: (3) or (18)
  integer,allocatable :: ngfft3_v1(:,:)
   ! ngfft3_v1(3, numv1)
   ! The FFT mesh used for each v1 potential (the one used to store data in the file).

 real(dp) :: acell(3),rprim(3,3),gprim(3,3)
   ! Used for slow FT.
   ! TODO: To be removed.

  real(dp),allocatable :: qpts(:,:)
   ! qpts(3,nqpt)
   ! List of q-points in reduced coordinates.

  ! Stuff needed for pawrhoij1
  !integer :: ntypat

  !integer,allocatable :: nlmn_type(:)
  !  nlmn_type(ntypat)= Number of (l,m,n) elements for the paw basis for each type of atom. Only used for reading.

  ! Fourier interpolation
  integer :: nrpt=0
  ! Number of real space points used for Fourier interpolation
  ! Note that the points are MPI distributed

  real(dp),allocatable :: rpt(:,:)
  ! rpt(3,nrpt)
  ! Real space points in canonical type coordinates.

  real(dp),allocatable :: v1scf_rpt(:,:,:,:,:)
  ! Wannier representation (real values)
  ! v1scf_rpt(nrpt, nfft , nspden, 3*natom)

  type(crystal_t) :: cryst
  ! Crystalline structure read from the the DVDV

  type(mpi_type) :: mpi_enreg
  ! Interl object used to call fourdp

 end type dvdb_t

 public :: dvdb_init              ! Initialize the object.
 public :: dvdb_open_read         ! Open the file in read-only mode.
 public :: dvdb_free              ! Close the file and release the memory allocated.
 public :: dvdb_print             ! Release memory.
 public :: dvdb_findq             ! Returns the index of the q-point.
 public :: dvdb_read_onev1        ! Read and return the DFPT potential for given (idir, ipert, iqpt).
 public :: dvdb_readsym_allv1     ! Read and return all the 3*natom DFPT potentials (either from file or symmetrized)
 public :: dvdb_list_perts        ! Check if all the (phonon) perts are available taking into account symmetries.
 public :: dvdb_check_v1sym       ! Check symmetries of the DFPT potentials.
 public :: dvdb_test_symmetries   ! Debugging tool used to test the symmetrization of the DFPT potentials.
 public :: dvdb_ftinterp_setup    ! Prepare the internal tables for Fourier interpolation.
 public :: dvdb_ftinterp_qpt      ! Fourier interpolation of potentials for given q-point
 public :: dvdb_merge_files       ! Merge a list of POT1 files.
!!***

CONTAINS
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_init
!! NAME
!!  dvdb_init
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
!!      eph,m_dvdb,mrgdv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_init(db, path, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_init'
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: path
 integer,intent(in) :: comm
 type(dvdb_t),intent(inout) :: db

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,timrev2=2
 integer :: iv1,ii,ierr,unt,fform,nqpt,iq,iq_found,cplex
 integer :: idir,ipert,my_rank
 character(len=500) :: msg
 type(hdr_type) :: hdr1,hdr_ref
!arrays
 integer,allocatable :: tmp_pos(:,:,:)
 real(dp),allocatable :: tmp_qpts(:,:)

!************************************************************************

 my_rank = xmpi_comm_rank(comm)
 db%path = path; db%comm = comm; db%iomode = IO_MODE_FORTRAN

 ! Master reads the header and builds useful tables
 if (my_rank == master) then

   if (open_file(path, msg, newunit=unt, form="unformatted", status="old", action="read") /= 0) then
     MSG_ERROR(msg)
   end if
   read(unt) db%version
   read(unt) db%numv1

   ! Get important dimensions from the first header and rewind the file.
   call hdr_fort_read(hdr_ref, unt, fform)
   ABI_CHECK(fform /= 0, sjoin("fform=0 while reading:", path))
   if (db%debug) call hdr_echo(hdr_ref,fform,4,unit=std_out)
   rewind(unt)
   read(unt)
   read(unt)

   ! The code below must be executed by the other procs if MPI.
   db%natom = hdr_ref%natom
   db%natom3 = 3 * hdr_ref%natom
   db%nspden = hdr_ref%nspden
   db%nsppol = hdr_ref%nsppol
   db%nspinor = hdr_ref%nspinor
   db%usepaw = hdr_ref%usepaw
   ABI_CHECK(db%usepaw == 0, "PAW not yet supported")

   ! TODO: Write function that returns mpert from natom!
   db%mpert = db%natom + 6

   ABI_MALLOC(tmp_qpts, (3, db%numv1))
   ABI_MALLOC(tmp_pos, (3, db%mpert, db%numv1))
   tmp_pos = 0

   ABI_MALLOC(db%cplex_v1, (db%numv1))
   ABI_MALLOC(db%ngfft3_v1, (3, db%numv1))

   nqpt = 0
   do iv1=1,db%numv1
     call hdr_fort_read(hdr1, unt, fform)
     ABI_CHECK(fform /= 0, sjoin("fform=0 while reading header of v1 potential of index:", itoa(iv1)))

     ! Save cplex and FFT mesh
     cplex = 2; if (hdr1%qptn(1)**2+hdr1%qptn(2)**2+hdr1%qptn(3)**2<1.d-14) cplex = 1
     db%cplex_v1(iv1) = cplex
     db%ngfft3_v1(:, iv1) = hdr1%ngfft(:3)

     ! Skip the records with v1.
     do ii=1,hdr1%nspden
       read(unt)
     end do

     ! check whether this q-point is already in the list.
     iq_found = 0
     do iq=1,nqpt
       if (all(abs(hdr1%qptn - tmp_qpts(:,iq)) < tol14)) then
         iq_found = iq; exit
       end if
     end do

     ! pertcase = idir + (ipert-1)*3 where ipert=iatom in the interesting cases
     idir = mod(hdr1%pertcase-1, 3) + 1
     ipert = (hdr1%pertcase - idir) / 3 + 1

     if (iq_found == 0) then
       nqpt = nqpt + 1
       tmp_qpts(:, nqpt) = hdr1%qptn
       iq_found = nqpt
     end if
     tmp_pos(idir, ipert, iq_found) = iv1

     call hdr_free(hdr1)
   end do

   ! Allocate arrays with correct shape.
   db%nqpt = nqpt
   ABI_MALLOC(db%qpts, (3, nqpt))
   db%qpts = tmp_qpts(:,1:nqpt)
   ABI_FREE(tmp_qpts)

   ABI_MALLOC(db%pos_dpq, (3, db%mpert, nqpt))
   db%pos_dpq = tmp_pos(:, :, 1:nqpt)
   ABI_FREE(tmp_pos)

   close(unt)
 end if

 ! Master broadcasts data.
 if (xmpi_comm_size(comm) > 1) then
   call xmpi_bcast(db%version, master, comm, ierr)
   call xmpi_bcast(db%numv1, master, comm, ierr)
   call xmpi_bcast(db%nqpt, master, comm, ierr)
   call hdr_bcast(hdr_ref, master, my_rank, comm)
   db%natom = hdr_ref%natom
   db%natom3 = 3 * hdr_ref%natom
   db%nspden = hdr_ref%nspden
   db%nsppol = hdr_ref%nsppol
   db%nspinor = hdr_ref%nspinor
   db%usepaw = hdr_ref%usepaw
   db%mpert = db%natom + 6

   if (my_rank /= master) then
     ABI_MALLOC(db%cplex_v1, (db%numv1))
     ABI_MALLOC(db%ngfft3_v1, (3, db%numv1))
     ABI_MALLOC(db%qpts, (3, db%nqpt))
     ABI_MALLOC(db%pos_dpq, (3, db%mpert, db%nqpt))
   end if

   call xmpi_bcast(db%cplex_v1, master, comm, ierr)
   call xmpi_bcast(db%ngfft3_v1, master, comm, ierr)
   call xmpi_bcast(db%qpts, master, comm, ierr)
   call xmpi_bcast(db%pos_dpq, master, comm, ierr)
 end if

 ! Init crystal_t from the hdr read from file.
 call crystal_from_hdr(db%cryst,hdr_ref,timrev2)
 call hdr_free(hdr_ref)

 ! Compute rprim, and gprimd. Used for slow FFT q--r if multiple shifts
 call mkradim(db%acell,db%rprim,db%cryst%rprimd)
 call matr3inv(db%rprim,db%gprim)

 ! Internal MPI_type needed for calling fourdp!
 call initmpi_seq(db%mpi_enreg)

end subroutine dvdb_init
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_open_read
!! NAME
!!  dvdb_open_read
!!
!! FUNCTION
!!  Open the file in read-only mode.
!!
!! INPUTS
!!   ngfft(18)=Info on the FFT used for the potentials. Note that ngfft
!!     is the mesh used by the parent and can differ from the one found in the file.
!!   comm=MPI communicator
!!
!! PARENTS
!!      m_phgamma,mrgdv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_open_read(db, ngfft, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_open_read'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 type(dvdb_t),intent(inout) :: db
!arrays
 integer,intent(in) :: ngfft(18)

!Local variables-------------------------------
!scalars
 integer :: nprocs
 character(len=500) :: msg

!************************************************************************

 if (db%rw_mode /= DVDB_NOMODE) then
   MSG_ERROR("DVDB should be in DVDB_NOMODE when open_read is called.")
 end if
 db%rw_mode = DVDB_READMODE

 nprocs = xmpi_comm_size(comm)

 ! Initialize tables to call fourdp in sequential
 db%ngfft = ngfft
 call init_distribfft_seq(db%mpi_enreg%distribfft,'c',ngfft(2),ngfft(3),'all')
 call init_distribfft_seq(db%mpi_enreg%distribfft,'f',ngfft(2),ngfft(3),'all')

 ! Open the file.
 select case (db%iomode)
 case (IO_MODE_FORTRAN)
   if (open_file(db%path, msg, newunit=db%fh, form="unformatted", status="old", action="read") /= 0) then
     MSG_ERROR(msg)
   end if
   read(db%fh)
   read(db%fh)
   db%current_fpos = 1

 case (IO_MODE_MPI)
   MSG_ERROR("MPI not coded")

 case default
   MSG_ERROR(sjoin("Unsupported iomode:", itoa(db%iomode)))
 end select

end subroutine dvdb_open_read
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_free
!! NAME
!!  dvdb_free
!!
!! FUNCTION
!! Close the file and release the memory allocated.
!!
!! PARENTS
!!      eph,m_dvdb,mrgdv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_free(db)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(dvdb_t),intent(inout) :: db

!************************************************************************

 ! integer arrays
 if (allocated(db%pos_dpq)) then
   ABI_FREE(db%pos_dpq)
 end if
 if (allocated(db%cplex_v1)) then
   ABI_FREE(db%cplex_v1)
 end if
 if (allocated(db%ngfft3_v1)) then
   ABI_FREE(db%ngfft3_v1)
 end if

 ! real arrays
 if (allocated(db%qpts)) then
   ABI_FREE(db%qpts)
 end if
 if (allocated(db%rpt)) then
   ABI_FREE(db%rpt)
 end if
 if (allocated(db%v1scf_rpt)) then
   ABI_FREE(db%v1scf_rpt)
 end if

 ! types
 call crystal_free(db%cryst)
 call destroy_mpi_enreg(db%mpi_enreg)

 ! Close the file but only if we have performed IO.
 if (db%rw_mode == DVDB_NOMODE) return

 select case (db%iomode)
 case (IO_MODE_FORTRAN)
   close(db%fh)
 case default
   MSG_ERROR(sjoin("Unsupported iomode:", itoa(db%iomode)))
 end select

end subroutine dvdb_free
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_print
!! NAME
!!  dvdb_print
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
!!      eph,m_dvdb,mrgdv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_print(db, header, unit, prtvol, mode_paral)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_print'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol,unit
 character(len=4),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header
 type(dvdb_t),intent(in) :: db

!Local variables-------------------------------
!scalars
 integer :: my_unt,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_unt =std_out; if (PRESENT(unit)) my_unt   =unit
 my_prtvol=0    ; if (PRESENT(prtvol)) my_prtvol=prtvol
 my_mode='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=' ==== Info on the dvdb% object ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 write(std_out,*)"Number of q-points: ",db%nqpt

end subroutine dvdb_print
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_get_pinfo
!! NAME
!!  dvdb_get_pinfo
!!
!! FUNCTION
!!  Return information on the perturbations available for a give q-point index.
!!
!! INPUTS
!!  iqpt=Index of the q-point
!!
!! OUTPUT
!!  nperts=Number of perturbations found.
!!  cplex=2 if potentials are complex, 1 for real
!!  pinfo(:,:)=Array with info on the perturbations present on file
!!     pinfo(1, ip) gives the `idir` index of the ip-th perturbation.
!!     pinfo(2, ip) gives the `ipert` index of the ip-th perturbation.
!!     pinfo(3, ip) gives `pertcase`=idir + (ipert-1)*3
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function dvdb_get_pinfo(db, iqpt, cplex, pinfo) result(nperts)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_get_pinfo'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(dvdb_t),intent(in) :: db
 integer,intent(in) :: iqpt
 integer,intent(out) :: cplex
!arrays
 integer,intent(out) :: pinfo(3,3*db%mpert)

!Local variables-------------------------------
!scalars
 integer :: idir,ipert,iv1

! *************************************************************************

 ! Get the number of perturbations computed for this iqpt
 nperts = 0; cplex = 0
 do ipert=1,db%natom ! natom selects atomic perturbations only.
    do idir=1,3
      iv1 = db%pos_dpq(idir,ipert,iqpt)
      if (iv1 /= 0) then
        nperts = nperts + 1
        pinfo(:, nperts) = [idir, ipert, idir + (ipert-1)*3]
        if (cplex == 0) cplex = db%cplex_v1(iv1)
        ABI_CHECK(cplex == db%cplex_v1(iv1), "cplex should be constant for given q!")
      end if
    end do
 end do

end function dvdb_get_pinfo
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_read_onev1
!! NAME
!!  dvdb_read_onev1
!!
!! FUNCTION
!!  Read the DFPT potential for the specified (idir, ipert, iqpt).
!!  Note that iqpt is the index in dvdb%qpts. Use dvdb_findq to
!!  get the index from the q-point in reduced coordinates.
!!
!! INPUTS
!!  idir=Direction of the perturbation
!!  ipert=Perturbation type.
!!  iqpt=Index of the q-point in dvdb%qpts
!!  cplex=1 if real, 2 if complex potentials.
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!
!! OUTPUT
!!  ierr=Non-zero if error.
!!  v1scf(cplex*nfft, nspden)=DFT potential associated to (idir, ipert, iqpt).
!!  msg=String with error message if ierr /= 0.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function dvdb_read_onev1(db, idir, ipert, iqpt, cplex, nfft, ngfft, v1scf, msg) result(ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_read_onev1'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idir,ipert,iqpt,cplex,nfft
 character(len=*),intent(out) :: msg
 type(dvdb_t),intent(inout) :: db
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(out) :: v1scf(cplex*nfft,db%nspden)

!Local variables-------------------------------
!scalars
 integer,parameter :: paral_kgb0=0
 integer :: iv1,ispden,nfftot_file,nfftot_out,ifft
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer :: ngfft_in(18),ngfft_out(18)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:),fftn3_distrib(:),ffti3_local(:)
 real(dp),allocatable :: v1r_file(:,:),v1g_in(:,:),v1g_out(:,:)

! *************************************************************************

 ! Consistency checks
 ierr = 1
 iv1 = db%pos_dpq(idir,ipert,iqpt)

 if (iv1 == 0) then
   write(msg,"(3(a,i0))")"Cannot find idir=",idir,", ipert=",ipert,", iqpt=",iqpt
   return
 end if

 if (cplex /= db%cplex_v1(iv1)) then
   write(msg,"(2(a,i0))")"Wrong cplex. Expecting=",db%cplex_v1(iv1),", received= ",cplex
   return
 end if

 call dvdb_seek(db, idir, ipert, iqpt)

 ierr = my_hdr_skip(db%fh, idir, ipert, db%qpts(:,iqpt))
 if (ierr /= 0)  then
   msg = sjoin("hdr_skip returned ierr:", itoa(ierr)); return
 end if

 ! Read v1 from file.
 nfftot_out = product(ngfft(:3))
 nfftot_file = product(db%ngfft3_v1(:3, iv1))

 if (all(ngfft(:3) == db%ngfft3_v1(:3, iv1))) then
   do ispden=1,db%nspden
     read(db%fh) (v1scf(ifft, ispden), ifft=1,cplex*nfftot_file)
   end do
 else
   ! Fourier interpolation
   !if (db%debug)
   MSG_ERROR("FFT interpolation of DFPT potentials must be tested.")
   ABI_MALLOC(v1r_file, (cplex*nfftot_file, db%nspden))
   do ispden=1,db%nspden
     read(db%fh) (v1r_file(ifft, ispden), ifft=1,cplex*nfftot_file)
   end do

   ! Call fourier_interpol to get v1scf on ngfft mesh.
   ngfft_in = ngfft; ngfft_out = ngfft
   ngfft_in(1:3) = db%ngfft3_v1(1:3, iv1); ngfft_out(1:3) = ngfft(1:3)
   ngfft_in(4:6) = ngfft_in(1:3); ngfft_out(4:6) = ngfft_out(1:3)
   ngfft_in(9:18) = 0; ngfft_out(9:18) = 0
   ngfft_in(10) = 1; ngfft_out(10) = 1

   call initmpi_seq(MPI_enreg_seq)
   ! Which one is coarse? Note that this part is not very robust and can fail!
   if (ngfft_in(2) * ngfft_in(3) < ngfft_out(2) * ngfft_out(3)) then
     call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfft_in(2),ngfft_in(3),'all')
     call init_distribfft_seq(MPI_enreg_seq%distribfft,'f',ngfft_out(2),ngfft_out(3),'all')
   else
     call init_distribfft_seq(MPI_enreg_seq%distribfft,'f',ngfft_in(2),ngfft_in(3),'all')
     call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfft_out(2),ngfft_out(3),'all')
   end if

   ABI_MALLOC(v1g_in,  (2, nfftot_file))
   ABI_MALLOC(v1g_out, (2, nfftot_out))

   call fourier_interpol(cplex,db%nspden,0,0,nfftot_file,ngfft_in,nfftot_out,ngfft_out,&
     paral_kgb0,MPI_enreg_seq,v1r_file,v1scf,v1g_in,v1g_out)

   ABI_FREE(v1g_in)
   ABI_FREE(v1g_out)
   ABI_FREE(v1r_file)
   call destroy_mpi_enreg(MPI_enreg_seq)
 end if

end function dvdb_read_onev1
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_readsym_allv1
!! NAME
!!  dvdb_readsym_allv1
!!
!! FUNCTION
!!  Read all 3*natom DFPT potentials for the given iqpt (only atomic perturbations).
!!
!!  The routine will:
!!     1) Reconstruct the potentials by symmetry if the DVDB contains less than 3*natom potentials.
!!     2) interpolate the data if the input FFT mesh defined by `ngfft` differs
!!        from the one used to store data in the file.
!!
!!  Note that iqpt is the index in dvdb%qpts. Use dvdb_findq to
!!  get the index from the q-point in reduced coordinates.
!!
!! INPUTS
!!  iqpt=Index of the q-point in dvdb%qpts
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  comm=MPI communicator
!!
!! OUTPUT
!!  cplex=1 if real, 2 if complex.
!!  v1scf(cplex, nfft, nspden, 3*natom)= v1scf potentials on the real-space FFT mesh for the 3*natom perturbations.
!!
!! PARENTS
!!      m_dvdb,m_phgamma,mrgdv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_readsym_allv1(db, iqpt, cplex, nfft, ngfft, v1scf, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_readsym_allv1'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqpt,nfft,comm
 integer,intent(out) :: cplex
 type(dvdb_t),target,intent(inout) :: db
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),allocatable,intent(out) :: v1scf(:,:,:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: ipc,npc,idir,ipert,pcase,my_rank,nproc,ierr
 character(len=500) :: msg
!arrays
 integer :: pinfo(3,3*db%mpert),pflag(3, db%natom)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 ! Get number of perturbations computed for this iqpt as well as cplex.
 npc = dvdb_get_pinfo(db, iqpt, cplex, pinfo)
 ABI_CHECK(npc /= 0, "npc == 0!")

 ABI_STAT_MALLOC(v1scf, (cplex,nfft,db%nspden,3*db%natom), ierr)
 ABI_CHECK(ierr==0, "OOM in v1scf")

 ! Master read all available perturbations and broadcasts data.
 if (my_rank == master) then
   do ipc=1,npc
     idir = pinfo(1,ipc); ipert = pinfo(2,ipc); pcase = pinfo(3, ipc)
     if (dvdb_read_onev1(db, idir, ipert, iqpt, cplex, nfft, ngfft, v1scf(:,:,:,pcase), msg) /= 0) then
       MSG_ERROR(msg)
     end if
   end do
 end if
 call xmpi_bcast(v1scf, master, comm, ierr)
 if (npc == 3*db%natom) return

 ! Perturbation are missing and we have to reconstruct them by symmetry.
 ! This is the common case when DFPT calculations are done for independent perturbations only.
 write(std_out,*)sjoin("Will use symmetries to recostruct:", itoa(3*db%natom - npc), "perturbations")

 ! 0 if pert is not available.
 ! 1 if pert is on file.
 ! 2 if pert has been reconstructed by symmetry.
 pflag = 0
 do ipc=1,npc
   pflag(pinfo(1,ipc), pinfo(2,ipc)) = 1
 end do

 call v1phq_complete(db%cryst,db%qpts(:,iqpt),ngfft,cplex,nfft,db%nspden,db%nsppol,db%mpi_enreg,db%symv1,pflag,v1scf)

end subroutine dvdb_readsym_allv1
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/v1phq_complete
!! NAME
!! v1phq_complete
!!
!! FUNCTION
!!  Use the symmetries of the little group of the q-point to reconstruct
!!  the first order potentials starting from an initial irreducible set.
!!
!! INPUTS
!!  cryst<crystal_t>=crystal structure parameters
!!  qpt(3)=q-point in reduced coordinates.
!!  ngfft(18)=Info of FFT grid.
!!  cplex=1 if real potentials (qpt==gamma), 2 if complex
!!  nfft=(effective) number of FFT grid points (for this proc).
!!  nspden=number of spin-density components
!!  nsppol=Number of independent spin polarizations
!!  mpi_enreg=information about MPI parallelization
!!  symv1=If True, the new potentials are symmetrized using the set of symmetries that leaves the
!!    perturbation invariant.
!!
!! SIDE EFFECTS
!!  pflag(3,natom)= For each atomic perturbation:
!!     0 if pert is not available. 1 if pert is available. 2 if pert has been reconstructed by symmetry.
!!     Initialized by the caller. Changed in output.
!!  v1scf(cplex*nfft,nspden,3*cryst%natom)=Array with first order potentials.
!!    in input: filled with the irreducible potentials (corresponding pflag set to 1)
!!    output: Contains full set of perturbations.
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine v1phq_complete(cryst,qpt,ngfft,cplex,nfft,nspden,nsppol,mpi_enreg,symv1,pflag,v1scf)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'v1phq_complete'
 use interfaces_32_util
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nspden,nsppol
 logical,intent(in) :: symv1
 type(crystal_t),intent(in) :: cryst
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(inout) :: pflag(3, cryst%natom)
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(inout) :: v1scf(cplex*nfft,nspden,3*cryst%natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: syuse0=0,rfmeth2=2,iout0=0,tim_fourdp0=0,iscf1=1
 integer :: idir,ipert,tsign,isym_eq,itirev_eq,ipert_eq !,itirev
 integer :: pcase,trev_q,idir_eq,pcase_eq,ispden,cnt!,trial
 integer :: i1,i2,i3,id1,id2,id3,n1,n2,n3,ind1,ind2,j1,j2,j3,l1,l2,l3,k1,k2,k3,nfftot
 real(dp) :: arg
 logical :: has_phase
 logical,parameter :: debug=.False.
 character(len=500) :: msg
!arrays
 integer :: symrel_eq(3,3),symrec_eq(3,3),symm(3,3),g0_qpt(3),l0(3),tsm1g(3)
 integer :: symq(4,2,cryst%nsym)
 real(dp) :: phnon1(2),tnon(3)
 real(dp),allocatable :: workg(:,:), workg_eq(:,:),v1g(:,:,:)

! *************************************************************************

 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3); nfftot = product(ngfft(1:3))
 ABI_CHECK(nfftot == nfft, "FFT parallelism not supported")
 id1 = n1/2+2; id2 = n2/2+2; id3 = n3/2+2

 ABI_MALLOC(v1g, (2,nfft,nspden))
 ABI_MALLOC(workg_eq, (2, nfft))
 ABI_MALLOC(workg, (2, nfft))

 ! Examine the symmetries of the q wavevector
 call littlegroup_q(cryst%nsym,qpt,symq,cryst%symrec,cryst%symafm,trev_q,prtvol=0)

pcase_loop: &
 do pcase=1,3*cryst%natom
   idir = mod(pcase-1, 3) + 1; ipert = (pcase - idir) / 3 + 1
   if (pflag(idir, ipert) /= 0) cycle ! This pcase is available

   ! Find symmetry which links to the perturbation requested (pcase)
   call find_symeq(cryst, idir, ipert, symq, pflag, ipert_eq, isym_eq, itirev_eq, g0_qpt)
   if (isym_eq == -1) then
     if (debug) write(std_out,*)"Cannot find isym eq for idir, ipert:", idir,ipert
     cycle pcase_loop
   end if

   ! set flag since we will reconstruct pcase from isym_eq.
   pflag(idir, ipert) = 2

   symrel_eq = cryst%symrel(:,:,isym_eq)
   symrec_eq = cryst%symrec(:,:,isym_eq)
   tsign = 3-2*itirev_eq

   ! Phase due to L0 + R^{-1}tau
   l0 = cryst%indsym(1:3,isym_eq,ipert_eq)
   call mati3inv(symrel_eq, symm); symm = transpose(symm)

   tnon = l0 + matmul(transpose(symrec_eq), cryst%tnons(:,isym_eq))
   has_phase = any(abs(tnon) > tol12)
   ! FIXME
   !ABI_CHECK(.not. has_phase, "has phase must be tested")
   if (has_phase) MSG_WARNING("has phase must be tested")

   workg = zero

   ! Reconstruct DFPT potential. Final results stored in v1g.
   !write(std_out,*)"Reconstructing idir, ipert",idir,ipert
   v1g = zero; cnt = 0
   do idir_eq=1,3
     if (symrec_eq(idir, idir_eq) == 0) cycle
     cnt = cnt + 1
     pcase_eq = idir_eq + (ipert_eq-1)*3
     if (debug) write(std_out,*)"idir_eq, ipert_eq, tsign",idir_eq, ipert_eq, tsign

     do ispden=1,nspden
       ! Get symmetric perturbation in G-space in workg_eq array.
       call fourdp(cplex,workg_eq,v1scf(:,ispden,pcase_eq),-1,mpi_enreg,nfft,ngfft,mpi_enreg%paral_kgb,tim_fourdp0)

       !call rotate_fqg(itirev_eq,symrec_eq,qpt,tnon,ngfft,nfft,nspden,workg_eq,workg)
       ind1=0
       do i3=1,n3
         do i2=1,n2
           do i1=1,n1
             ind1=ind1+1

             ! Get location of G vector (grid point) centered at 0 0 0
             l1 = i1-(i1/id1)*n1-1
             l2 = i2-(i2/id2)*n2-1
             l3 = i3-(i3/id3)*n3-1

             ! Get rotated G vector Gj for each symmetry element
             ! -- here we use the TRANSPOSE of symrel_eq; assuming symrel_eq expresses
             ! the rotation in real space, the transpose is then appropriate
             ! for G space symmetrization (p. 1172d,e of notes, 2 June 1995).
             j1 = tsign * (symrel_eq(1,1)*l1+symrel_eq(2,1)*l2+symrel_eq(3,1)*l3)
             j2 = tsign * (symrel_eq(1,2)*l1+symrel_eq(2,2)*l2+symrel_eq(3,2)*l3)
             j3 = tsign * (symrel_eq(1,3)*l1+symrel_eq(2,3)*l2+symrel_eq(3,3)*l3)

             ! FIXME :TO BE CLARIFIED:
             ! We are not working on the G-sphere thus SG may be outside
             ! of the box. This check is not done in irrzg!!!
             if ( (j1 > n1/2 .or. j1 < -(n1-1)/2) .or. &
                  (j2 > n2/2 .or. j1 < -(n2-1)/2) .or. &
                  (j3 > n3/2 .or. j3 < -(n3-1)/2) ) then
               !write(std_out,*)"got it"
               workg(:, ind1) = zero; cycle
             end if

             tsm1g = [j1,j2,j3] ! +- S^{-1} G

             ! Map into [0,n-1] and then add 1 for array index in [1,n]
             k1=1+mod(n1+mod(j1,n1),n1)
             k2=1+mod(n2+mod(j2,n2),n2)
             k3=1+mod(n3+mod(j3,n3),n3)

             ! Get linear index of rotated point Gj
             ind2 = k1+n1*((k2-1)+n2*(k3-1))

             if (has_phase) then
               ! compute exp(-2*Pi*I*G dot tau) using original G
               ! NB: this phase is same as that in irrzg and phnons1, and corresponds
               ! to complex conjugate of phase from G to Gj;
               ! we use it immediately below, to go _to_ workg_eq(ind1)
               arg = two_pi * dot_product(qpt + tsm1g, tnon)
               phnon1(1) = cos(arg); phnon1(2) = -sin(arg)

               ! rho(Strans*G)=exp(2*Pi*I*(G) dot tau_S) rho(G)
               workg(1, ind1) = phnon1(1) * workg_eq(1, ind2) - phnon1(2) * workg_eq(2, ind2)
               workg(2, ind1) = phnon1(1) * workg_eq(2, ind2) + phnon1(2) * workg_eq(1, ind2)
             else
               workg(1, ind1) = workg_eq(1, ind2)
               workg(2, ind1) = workg_eq(2, ind2)
             end if

             ! Take complex conjugate if time-reversal is used.
             if (tsign == -1) workg(2, ind1) = -workg(2, ind1)
           end do
         end do
       end do

       v1g(:,:,ispden) = v1g(:,:,ispden) + workg * symrec_eq(idir, idir_eq)
     end do ! ispden
   end do ! idir_eq
   if (debug) write(std_out,*)"Used ",cnt," equivalent perturbations"

   ! Get potential in real space (results in v1scf)
   do ispden=1,nspden
     call fourdp(cplex,v1g(:,:,ispden),v1scf(:,ispden,pcase),+1,mpi_enreg,nfft,ngfft,mpi_enreg%paral_kgb,tim_fourdp0)

     ! IS(q) = q + G0
     ! we want q so we have to multiply by exp(iG0r) in real space.
     if (any(g0_qpt /= 0)) then
       ABI_CHECK(cplex==2, "cplex == 1")
       if (debug) write(std_out,*)"Found not zero g0_qpt for idir, ipert:",idir,ipert
       call times_eigr(g0_qpt, ngfft, nfft, 1, v1scf(:,ispden,pcase))
     end if
   end do

   if (symv1) then
     if (debug) write(std_out,*)"Calling v1phq_symmetrize"
     call v1phq_symmetrize(cryst,idir,ipert,symq,ngfft,cplex,nfft,nspden,nsppol,mpi_enreg,v1scf(:,:,pcase))
   end if
 end do pcase_loop

 ABI_FREE(v1g)
 ABI_FREE(workg)
 ABI_FREE(workg_eq)

 ! Handle possible error.
 if (any(pflag == 0)) then
   write(std_out,"(2a)")&
     "The following perturbations cannot be recostructed by symmetry for q-point: ",trim(ktoa(qpt))
   do ipert=1,cryst%natom
     do idir=1,3
        if (pflag(idir, ipert) == 0) write(std_out,"(2(a,i0))")"idir= ",idir,", ipert= ",ipert
     end do
   end do
   write(msg,"(5a)")&
     "Cannot recostruct all 3*natom atomic perturbations from file",ch10,&
     "This usually happens when the DVDB does not contain all the independent perturbations for this q-point",ch10,&
     "See message above for further info"
   MSG_ERROR(msg)
 end if

end subroutine v1phq_complete
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/find_symeq
!! NAME
!! find_symeq
!!
!! FUNCTION
!!  Find symmetry which links to the perturbation specified by (idir, ipert)
!!
!! INPUTS
!!  cryst<crystal_t>=crystal structure parameters
!!  idir=Direction of the perturbation
!!  ipert=Perturbation type.
!!  symq(4,2,nsym)=Table produced by littlegroup_q
!!  pflag(3,natom)= For each atomic perturbation:
!!     0 if pert is not available. 1 if pert is available. 2 if pert can been reconstructed by symmetry.
!!
!! OUTPUT
!!  ipert_eq
!!  isym_eq
!!  itirev_eq
!!  g0_qpt(3)
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine find_symeq(cryst, idir, ipert, symq, pflag, ipert_eq, isym_eq, itirev_eq, g0_qpt)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'find_symeq'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idir,ipert
 integer,intent(out) :: ipert_eq,isym_eq,itirev_eq
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: symq(4,2,cryst%nsym),pflag(3,cryst%natom)
 integer,intent(out) :: g0_qpt(3)

!Local variables-------------------------------
!scalars
 integer :: isym,idir_eq,ip,itirev

! *************************************************************************

 isym_eq = -1; ipert_eq = -1

symloop: &
 do itirev=1,2
   itirev_eq = itirev
   do isym=1,cryst%nsym

     ! Check that isym preserves the q-point
     !if (symq(4,itirev,isym) /= 1 .or. any(symq(1:3,itirev,isym) /= 0)) cycle
     if (symq(4,itirev,isym) /= 1) cycle ! .or. any(symq(1:3,itirev,isym) /= 0)) cycle
     g0_qpt = symq(1:3,itirev,isym)

     do ip=1,cryst%natom
       !if (.not. cryst%indsym(4,isym,ip) == ipert) cycle
       if (.not. cryst%indsym(4,isym,ipert) == ip) cycle
       isym_eq = isym; ipert_eq = ip
       do idir_eq=1,3
         if (idir_eq == idir .and. ip == ipert .and. cryst%symrec(idir,idir_eq,isym) /= 0) isym_eq = -1
         if (cryst%symrec(idir,idir_eq,isym) /= 0 .and. pflag(idir_eq, ip) == 0) then
         !if (cryst%symrel(idir,idir_eq,isym) /= 0 .and. pflag(idir_eq, ip) == 0) then
           !if (idir_eq == idir .and. ip == ipert) cycle
           isym_eq = -1
         end if
       end do
       if (isym_eq /= -1) exit symloop
     end do
   end do
 end do symloop

 if (isym_eq == -1) then
   ipert_eq = -1; itirev_eq = -1
 end if

end subroutine find_symeq
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/v1phq_rotate
!! NAME
!! v1phq_rotate
!!
!! FUNCTION
!!  Reconstruct the DFPT potential for q in BZ starting from the value in the IBZ.
!!
!! INPUTS
!!  cryst<crystal_t>=crystal structure parameters
!!  qpt_ibz(3)=q-point in the IBZ in reduced coordinates.
!!  ngfft=array of dimensions for different FFT grids
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  nfft=(effective) number of FFT grid points (for this proc) for the "fine" grid (see NOTES in respfn.F90)
!!  nspden=number of spin-density components
!!  nsppol=Number of independent spin polarizations
!!  mpi_enreg=information about MPI parallelization
!!  v1r_qibz(cplex*nfft,nspden,3*cryst%natom)=Array with first order potentials in real space
!!   for the irreducible q-point `qpt_ibz`
!!
!! OUTPUT
!!  v1r_qbz(cplex*nfft,nspden,3*cryst%natom)=Array with first order potentials in real space for the q-point in the BZ
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine v1phq_rotate(cryst,qpt_ibz,isym,itimrev,g0q,ngfft,cplex,nfft,nspden,nsppol,mpi_enreg,v1r_qibz,v1r_qbz)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'v1phq_rotate'
 use interfaces_32_util
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isym,itimrev,cplex,nfft,nspden,nsppol
 type(crystal_t),intent(in) :: cryst
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: g0q(3),ngfft(18)
 real(dp),intent(in) :: qpt_ibz(3)
 real(dp),intent(inout) :: v1r_qibz(cplex*nfft,nspden,3*cryst%natom)
 real(dp),intent(out) :: v1r_qbz(cplex*nfft,nspden,3*cryst%natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourdp0=0
 integer :: natom3,mu,ispden,idir,ipert,idir_eq,ipert_eq,mu_eq,cnt,tsign
!arrays
 integer :: symrec_eq(3,3),sm1(3,3),l0(3) !g0_qpt(3), symrel_eq(3,3),
 real(dp) :: tnon(3)
 real(dp),allocatable :: v1g_qibz(:,:,:),workg(:,:),v1g_mu(:,:)

! *************************************************************************

 ABI_UNUSED(nsppol)
 ABI_CHECK(cplex == 2, "cplex != 2")

 natom3 = 3 * cryst%natom
 tsign = 3-2*itimrev

 ! Compute IBZ potentials in G-space (results in v1g_qibz)
 ABI_MALLOC(v1g_qibz, (2*nfft,nspden,natom3))
 do mu=1,natom3
   do ispden=1,nspden
     call fourdp(cplex,v1g_qibz(:,ispden,mu),v1r_qibz(:,ispden,mu),-1,mpi_enreg,nfft,ngfft,mpi_enreg%paral_kgb,tim_fourdp0)
   end do
 end do

 ABI_MALLOC(workg, (2*nfft,nspden))
 ABI_MALLOC(v1g_mu, (2*nfft,nspden))

 ! For each perturbation:
 ! FIXME: This is wrong if natom > 1
 symrec_eq = cryst%symrec(:,:,isym)
 call mati3inv(symrec_eq, sm1); sm1 = transpose(sm1)

 do mu=1,natom3
   idir = mod(mu-1, 3) + 1; ipert = (mu - idir) / 3 + 1

   ! Phase due to L0 + R^{-1}tau
   l0 = cryst%indsym(1:3,isym,ipert)
   tnon = l0 + matmul(transpose(symrec_eq), cryst%tnons(:,isym))
   ! FIXME
   !ABI_CHECK(all(abs(tnon) < tol12), "tnon!")
   if (.not.all(abs(tnon) < tol12)) MSG_WARNING("tnon!")

   ipert_eq = cryst%indsym(4,isym,ipert)

   v1g_mu = zero; cnt = 0
   do idir_eq=1,3
     if (symrec_eq(idir, idir_eq) == 0) cycle
     mu_eq = idir_eq + (ipert_eq-1)*3
     cnt = cnt + 1

     ! Rotate in G-space and accumulate in workg
     call rotate_fqg(itimrev,sm1,qpt_ibz,tnon,ngfft,nfft,nspden,v1g_qibz(:,:,mu_eq),workg)
     v1g_mu = v1g_mu + workg * symrec_eq(idir, idir_eq)
   end do ! idir_eq
   ABI_CHECK(cnt /= 0, "cnt should not be zero!")

   ! Transform to real space and take into account a possible shift.
   ! (results are stored in v1r_qbz)
   do ispden=1,nspden
     call fourdp(cplex,v1g_mu(:,ispden),v1r_qbz(:,ispden,mu),+1,mpi_enreg,nfft,ngfft,mpi_enreg%paral_kgb,tim_fourdp0)
     call times_eigr(-g0q, ngfft, nfft, 1, v1r_qbz(:,ispden,mu))
     !call times_eigr(tsign * g0q, ngfft, nfft, 1, v1r_qbz(:,ispden,mu))
   end do

   !call v1phq_symmetrize(cryst,idir,ipert,symq,ngfft,cplex,nfft,nspden,nsppol,mpi_enreg,v1r)
 end do ! mu

 ABI_FREE(workg)
 ABI_FREE(v1g_mu)
 ABI_FREE(v1g_qibz)

end subroutine v1phq_rotate
!!***

!!****f* m_dvdb/v1phq_symmetrize
!! NAME
!! v1phq_symmetrize
!!
!! FUNCTION
!!  Enforce spatial-symmetry on the DFPT potential.
!!
!! INPUTS
!! cryst<crystal_t>=crystal structure parameters
!! idir=Direction of the perturbation
!! ipert=Perturbation type.
!! symq(4,2,nsym)= Table computed by littlegroup_q.
!!   three first numbers define the G vector;
!!   fourth number is zero if the q-vector is not preserved, is 1 otherwise
!!   second index is one without time-reversal symmetry, two with time-reversal symmetry
!! ngfft=array of dimensions for different FFT grids
!! cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!! nfft=(effective) number of FFT grid points (for this proc) for the "fine" grid (see NOTES in respfn.F90)
!! nspden=number of spin-density components
!! mpi_enreg=information about MPI parallelization
!!
!! SIDE EFFECTS
!!  v1r(cplex*nfft,nspden,3*cryst%natom)=Array with first order potentials in real space
!!   Symmetrized in output.
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine v1phq_symmetrize(cryst,idir,ipert,symq,ngfft,cplex,nfft,nspden,nsppol,mpi_enreg,v1r)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'v1phq_symmetrize'
 use interfaces_41_geometry
 use interfaces_56_recipspace
 use interfaces_67_common
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idir,ipert,cplex,nfft,nspden,nsppol
 type(crystal_t),intent(in) :: cryst
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: symq(4,2,cryst%nsym),ngfft(18)
 real(dp),intent(inout) :: v1r(cplex*nfft,nspden)

!Local variables-------------------------------
 integer,parameter :: syuse0=0,rfmeth2=2,iscf1=1 !iout0=0,
 integer :: nsym1,nfftot
!arrays
 integer :: symafm1(cryst%nsym),symrel1(3,3,cryst%nsym),symrc1(3,3,cryst%nsym)
 integer,allocatable :: irrzon1(:,:,:),indsy1(:,:,:)
 real(dp) :: tnons1(3,cryst%nsym)
 real(dp),allocatable :: phnons1(:,:,:),v1g(:,:)

! *************************************************************************

 if (cryst%nsym == 1) return

 nfftot = product(ngfft(1:3))
 ABI_CHECK(nfft == nfftot, "MPI-FFT not coded")

 ! Symmetrize (copied from dfpt_looppert)
 ! Determines the set of symmetries that leaves the perturbation invariant.
 call littlegroup_pert(cryst%gprimd,idir,cryst%indsym,std_out,ipert,cryst%natom,cryst%nsym,nsym1,rfmeth2,&
   cryst%symafm,symafm1,symq,cryst%symrec,cryst%symrel,symrel1,syuse0,cryst%tnons,tnons1)

 ! Set up corresponding symmetry data
 ABI_MALLOC(irrzon1, (nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4)))
 ABI_MALLOC(phnons1, (2,nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4)))
 ABI_MALLOC(indsy1,(4,nsym1,cryst%natom))

 call setsym(indsy1,irrzon1,iscf1,cryst%natom,nfft,ngfft,nspden,nsppol,&
   nsym1,phnons1,symafm1,symrc1,symrel1,tnons1,cryst%typat,cryst%xred)

 !if (psps%usepaw==1) then
 !  ! Allocate/initialize only zarot in pawang1 datastructure
 !  call pawang_init(pawang1,0,pawang%l_max-1,0,nsym1,0,1,0,0,0)
 !  call setsymrhoij(gprimd,pawang1%l_max-1,pawang1%nsym,0,rprimd,symrc1,pawang1%zarot)
 !end if

 ! FIXME Be careful here because symrhg was written for densities!
 ABI_CHECK(nsppol == 1 .and. nspden == 1, "symrhg was written for densities, not for potentials")

 ABI_MALLOC(v1g, (2,nfft))
 call symrhg(cplex,cryst%gprimd,irrzon1,mpi_enreg,nfft,nfftot,ngfft,nspden,nsppol,nsym1,&
    mpi_enreg%paral_kgb,phnons1,v1g,v1r,cryst%rprimd,symafm1,symrel1)

 ABI_FREE(irrzon1)
 ABI_FREE(phnons1)
 ABI_FREE(indsy1)
 ABI_FREE(v1g)

end subroutine v1phq_symmetrize
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/rotate_fqg
!! NAME
!!  rotate_fqg
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine rotate_fqg(itirev,symm,qpt,tnon,ngfft,nfft,nspden,infg,outfg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rotate_fqg'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: itirev,nfft,nspden
!scalars
 integer,intent(in) :: ngfft(18)
!arrays
 integer,intent(in) :: symm(3,3)
 real(dp),intent(in) :: qpt(3),tnon(3)
 real(dp),intent(in) :: infg(2,nfft,nspden)
 real(dp),intent(out) :: outfg(2,nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,id1,id2,id3,n1,n2,n3,ind1,ind2,j1,j2,j3,l1,l2,l3,k1,k2,k3,nfftot,isp,tsign
 real(dp) :: arg
 logical :: has_phase
!arrays
 integer :: tsg(3)
 real(dp) :: phnon1(2)

! *************************************************************************

 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3); nfftot = product(ngfft(1:3))
 ABI_CHECK(nfftot == nfft, "FFT parallelism not supported")
 id1 = n1/2+2; id2 = n2/2+2; id3=n3/2+2

 ABI_CHECK(any(itirev == [1,2]), "Wrong itirev")
 tsign = 3-2*itirev; has_phase = any(abs(tnon) > tol12)

 do isp=1,nspden
   ind1 = 0
   do i3=1,n3
     do i2=1,n2
       do i1=1,n1
         ind1 = ind1 + 1

         ! Get location of G vector (grid point) centered at 0 0 0
         l1 = i1-(i1/id1)*n1-1
         l2 = i2-(i2/id2)*n2-1
         l3 = i3-(i3/id3)*n3-1

         ! Get rotated G vector. IS(G)
         j1 = tsign * (symm(1,1)*l1+symm(1,2)*l2+symm(1,3)*l3)
         j2 = tsign * (symm(2,1)*l1+symm(2,2)*l2+symm(2,3)*l3)
         j3 = tsign * (symm(3,1)*l1+symm(3,2)*l2+symm(3,3)*l3)

         ! FIXME :TO BE CLARIFIED:
         ! We are not working on the G-sphere thus SG may be outside
         ! of the box. This check is not done in irrzg!!!
         if ( (j1 > n1/2 .or. j1 < -(n1-1)/2) .or. &
              (j2 > n2/2 .or. j1 < -(n2-1)/2) .or. &
              (j3 > n3/2 .or. j3 < -(n3-1)/2) ) then
           !write(std_out,*)"got it"
           outfg(:,ind1,isp) = zero
           cycle
         end if

         tsg = [j1,j2,j3] ! +- S^{-1} G

         ! Map into [0,n-1] and then add 1 for array index in [1,n]
         k1=1+mod(n1+mod(j1,n1),n1)
         k2=1+mod(n2+mod(j2,n2),n2)
         k3=1+mod(n3+mod(j3,n3),n3)

         ! Get linear index of rotated point Gj
         ind2 = k1+n1*((k2-1)+n2*(k3-1))

         if (has_phase) then
           ! compute exp(-2*Pi*I*G dot tau) using original G
           arg = two_pi * dot_product(qpt + tsg, tnon)
           phnon1(1) = cos(arg); phnon1(2) =-sin(arg)

           ! rho(Strans*G)=exp(2*Pi*I*(G) dot tau_S) rho(G)
           outfg(1, ind1, isp) = phnon1(1) * infg(1, ind2, isp) - phnon1(2) * infg(2, ind2, isp)
           outfg(2, ind1, isp) = phnon1(1) * infg(2, ind2, isp) + phnon1(2) * infg(1, ind2, isp)
         else
           outfg(1, ind1, isp) = infg(1, ind2, isp)
           outfg(2, ind1, isp) = infg(2, ind2, isp)
         end if

         ! Take complex conjugate if time-reversal is used.
         if (tsign == -1) outfg(2, ind1, isp) = -outfg(2, ind1, isp)
       end do
     end do
   end do
 end do ! isp

end subroutine rotate_fqg
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_ftinterp_setup
!! NAME
!!  dvdb_ftinterp_setup
!!
!! FUNCTION
!!  Prepare the internal tables for the Fourier interpolation of the DFPT potentials.
!!
!! INPUTS
!!  ngqpt(3)=Divisions of the ab-initio q-mesh.
!!  nqshift=Number of shifts used to generated the ab-initio q-mesh.
!!  qshift(3,nqshift)=The shifts of the ab-initio q-mesh.
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  comm=MPI communicator
!!
!! PARENTS
!!      m_phgamma,mrgdv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_ftinterp_setup(db,ngqpt,nqshift,qshift,nfft,ngfft,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_ftinterp_setup'
 use interfaces_51_manage_mpi
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqshift,nfft,comm
 type(dvdb_t),target,intent(inout) :: db
!arrays
 integer,intent(in) :: ngqpt(3),ngfft(18)
 real(dp),intent(in) :: qshift(3,nqshift)

!Local variables-------------------------------
!scalars
 integer,parameter :: chksymbreak0=0,iout0=0,iscf2=2,brav1=1,sppoldbl1=1,timrev1=1,tim_fourdp0=0,parak_kgb0=0
 integer :: my_qptopt,iq_ibz,nqibz,iq_bz,nqbz,nqpt_computed
 integer :: ii,my_nqshift,iq_dvdb,cplex_qibz,ispden,mu,irpt
 integer :: iqst,nqst,itimrev,tsign,isym,ix,iy,iz,nq1,nq2,nq3,r1,r2,r3
 integer :: nproc,my_rank,ifft,cnt,ierr
 real(dp) :: qptrlen,dksqmax,phre,phim
 logical :: isirr_q
 type(crystal_t),pointer :: cryst
 type(mpi_type) :: mpi_enreg_seq
!arrays
 integer,parameter :: vacuum0(3)=[0,0,0]
 integer :: qptrlatt(3,3),g0q(3),ngfft_qspace(18)
 integer,allocatable :: indqq(:,:),iperm(:),bz2ibz_sort(:)
 real(dp) :: my_qshift(3,210),qpt_bz(3),shift(3) !,qpt_ibz(3)
 real(dp),allocatable :: qibz(:,:),qbz(:,:),wtq(:),emiqr(:,:)
 real(dp),allocatable :: v1r_qibz(:,:,:,:),v1r_qbz(:,:,:,:) !,all_v1qr(:,:,:,:)

! *************************************************************************

 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 cryst => db%cryst
 nq1= ngqpt(1); nq2 = ngqpt(2); nq3 = ngqpt(3)

 my_qptopt = 1 !; if (present(qptopt)) my_qptopt = qptopt

 ! Generate the q-mesh by finding the IBZ and the corresponding weights.
 qptrlatt = 0
 do ii=1,3
   qptrlatt(ii,ii) = ngqpt(ii)
 end do
 my_nqshift = nqshift ! Be careful as getkgrid expects shiftk(3,210).
 ABI_CHECK(my_nqshift>0 .and. my_nqshift<=210, "nqshift must be in [1,210]")

 ! First call to getkgrid to obtain nqibz.
 my_qshift=zero; my_qshift(:,1:nqshift) = qshift
 ABI_MALLOC(qibz, (3, 0))
 ABI_MALLOC(wtq, (0))
 call getkgrid(chksymbreak0,iout0,iscf2,qibz,my_qptopt,qptrlatt,qptrlen,cryst%nsym,0,nqibz,&
   my_nqshift,cryst%nsym,cryst%rprimd,my_qshift,cryst%symafm,cryst%symrel,vacuum0,wtq)
 ABI_FREE(qibz)
 ABI_FREE(wtq)

 ! Recall getkgrid to get qibz and wtq.
 ABI_MALLOC(qibz,(3,nqibz))
 ABI_MALLOC(wtq,(nqibz))

 call getkgrid(chksymbreak0,iout0,iscf2,qibz,my_qptopt,qptrlatt,qptrlen,cryst%nsym,nqibz,nqpt_computed,&
   my_nqshift,cryst%nsym,cryst%rprimd,my_qshift,cryst%symafm,cryst%symrel,vacuum0,wtq,fullbz=qbz)
 nqbz = size(qbz, dim=2)

 write(std_out,*)"Irreducible q-points:"
 do iq_ibz=1,nqibz
   write(std_out,*)trim(ktoa(qibz(:,iq_ibz))),wtq(iq_ibz)*nqbz
 end do

 ABI_CHECK(nqbz == product(ngqpt) * nqshift, "nqbz /= product(ngqpt) * nqshift")

 !nqbz = product(ngqpt)
 ABI_FREE(qbz)
 ABI_MALLOC(qbz, (3, nqbz))
 ii = 0
 do iz=0,nq3-1
   do iy=0,nq2-1
     do ix=0,nq1-1
       ii = ii + 1
       qbz(:, ii) = [ix/dble(nq1), iy/dble(nq2), iz/dble(nq3)]
       call wrap2_pmhalf([ix/dble(nq1), iy/dble(nq2), iz/dble(nq3)], qbz(:,ii), shift)
     end do
   end do
 end do

 ! Compute real-space points.
 ! Use the following indexing (N means ngfft of the adequate direction)
 ! 0 1 2 3 ... N/2    -(N-1)/2 ... -1    <= gc
 ! 1 2 3 4 ....N/2+1  N/2+2    ...  N    <= index ig
 db%nrpt = nqbz
 ABI_MALLOC(db%rpt, (3, db%nrpt))
 ii = 0
 do iz=1,nq3
   r3 = ig2gfft(iz,nq3) !* nq3
   do iy=1,nq2
     r2 = ig2gfft(iy,nq2) !* nq2
     do ix=1,nq1
       r1 = ig2gfft(ix,nq1) !* nq1
       ii = ii + 1
       db%rpt(:,ii) = [r1, r2, r3]
     end do
   end do
 end do

 ! Find correspondence BZ --> IBZ. Note:
 !   * time-reversal is always used for phonons
 !   * we use symrec instead of symrel
 ABI_MALLOC(indqq, (nqbz*sppoldbl1,6))

 call listkk(dksqmax,cryst%gmet,indqq,qibz,qbz,nqibz,nqbz,cryst%nsym,&
   sppoldbl1,cryst%symafm,cryst%symrec,timrev1,use_symrec=.True.)

 if (dksqmax > tol12) then
   MSG_BUG("Something wrong in the generation of the q-points in the BZ! Cannot map BZ --> IBZ")
 end if

 ! Construct sorted mapping BZ --> IBZ to speedup qbz search below.
 ABI_MALLOC(iperm, (nqbz))
 ABI_MALLOC(bz2ibz_sort, (nqbz))
 iperm = [(ii, ii=1,nqbz)]
 bz2ibz_sort = indqq(:,1)
 call sort_int(nqbz,bz2ibz_sort,iperm)

 ABI_MALLOC(emiqr, (2,db%nrpt))
 ABI_STAT_MALLOC(db%v1scf_rpt,(2,db%nrpt, nfft, db%nspden, db%natom3), ierr)
 ABI_CHECK(ierr==0, 'out of memory in db%v1scf_rpt')
 db%v1scf_rpt = zero

 ABI_MALLOC(v1r_qbz, (2, nfft, db%nspden, db%natom3))
 !v1r_qbz = huge(one)

 iqst = 0
 do iq_ibz=1,nqibz

   iq_dvdb = dvdb_findq(db, qibz(:,iq_ibz))
   if (iq_dvdb == -1) then
     MSG_ERROR(sjoin("Cannot find q-point:", ktoa(qibz(:,iq_ibz)), "in DVDB file"))
   end if
   !write(std_out,*)sjoin("qpt irred:",ktoa(qibz(:,iq_ibz)))

   ! Get potentials for this IBZ q-point on the real-space FFT mesh.
   ! This call allocates v1r_qibz(cplex_qibz, nfft, nspden, 3*natom)
   call dvdb_readsym_allv1(db, iq_dvdb, cplex_qibz, nfft, ngfft, v1r_qibz, comm)

   ! Find number of symmetric q-ponts associated to iq_ibz
   nqst = 0
   do ii=iqst+1,nqbz
     if (bz2ibz_sort(ii) /= iq_ibz) exit
     nqst = nqst + 1
   end do
   ABI_CHECK(nqst > 0 .and. bz2ibz_sort(iqst+1) == iq_ibz, "Wrong iqst")
   if (abs(nqst - wtq(iq_ibz) * nqbz) > tol12) then
     write(std_out,*)nqst, wtq(iq_ibz) * nqbz
     MSG_ERROR("Stop")
   end if

   ! Reconstruct by symmetry the potentials for the star of this q-point, perform FT and accumulate
   ! Be careful with the gamma point.
   do ii=1,nqst
     iqst = iqst + 1
     iq_bz = iperm(iqst)
     ABI_CHECK(iq_ibz == indqq(iq_bz,1), "iq_ibz !/ ind qq(1)")
     isym = indqq(iq_bz,2); itimrev = indqq(iq_bz,6) + 1; g0q = indqq(iq_bz,3:5) ! IS(q_ibz) + g0q = q_bz
     tsign = 3-2*itimrev

     qpt_bz = qbz(:, iq_bz)
     !write(std_out,*)"  treating:",trim(ktoa(qpt_bz))
     isirr_q = (isym == 1 .and. itimrev == 1 .and. all(g0q == 0))
     !ABI_CHECK(all(g0q == 0), "g0q /= 0")

     if (cplex_qibz == 1) then
       ! Gamma point.
       ABI_CHECK(nqst == 1, "cplex_qibz == 1 and nq nqst /= 1 (should be gamma)")
       ABI_CHECK(all(g0q == 0), "gamma point with g0q /= 0")

       ! Slow FT.
       cnt = 0
       do mu=1,db%natom3
         do ispden=1,db%nspden
           do irpt=1,db%nrpt
             cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
             do ifft=1,nfft
               db%v1scf_rpt(1, irpt, ifft, ispden, mu) = db%v1scf_rpt(1, irpt, ifft, ispden, mu) + &
                  v1r_qibz(1, ifft, ispden, mu)
             end do
           end do
         end do
       end do

     else
       ! Get the periodic part of the potential in BZ (v1r_qbz)
       if (isirr_q) then
         !write(std_out,*)sjoin("qpt irred:",ktoa(qpt_bz))
         v1r_qbz = v1r_qibz
       else
         call v1phq_rotate(cryst, qibz(:,iq_ibz), isym, itimrev, g0q,&
           ngfft, cplex_qibz, nfft, db%nspden, db%nsppol, db%mpi_enreg, v1r_qibz, v1r_qbz)
         !v1r_qbz = zero
         !v1r_qbz = v1r_qibz

         !call times_eigr(-tsign * g0q, ngfft, nfft, db%nspden*db%natom3, v1r_qbz)
         !call times_eigr(+tsign * g0q, ngfft, nfft, db%nspden*db%natom3, v1r_qbz)
         !if (itimrev == 2) v1r_qbz(2,:,:,:) = -v1r_qbz(2,:,:,:)
         !call times_eigr(-tsign * g0q, ngfft, nfft, db%nspden*db%natom3, v1r_qbz)
         !call times_eigr(+tsign * g0q, ngfft, nfft, db%nspden*db%natom3, v1r_qbz)
       end if

       ! Multiply by e^{iqpt_bz.r}
       call times_eikr(qpt_bz, ngfft, nfft, db%nspden*db%natom3, v1r_qbz)

       ! Compute FT phases for this qpt_bz.
       call calc_eiqr(-qpt_bz,db%nrpt,db%rpt,emiqr)

       ! Slow FT.
       cnt = 0
       do mu=1,db%natom3
         do ispden=1,db%nspden
           do irpt=1,db%nrpt
             cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
             phre = emiqr(1,irpt); phim = emiqr(2,irpt)

             do ifft=1,nfft
               db%v1scf_rpt(1,irpt, ifft, ispden, mu) = db%v1scf_rpt(1, irpt, ifft, ispden, mu) + &
                  phre * v1r_qbz(1, ifft, ispden, mu) - phim * v1r_qbz(2, ifft, ispden, mu)

               db%v1scf_rpt(2, irpt, ifft, ispden, mu) = db%v1scf_rpt(2, irpt, ifft, ispden, mu) + &
                  phre * v1r_qbz(2, ifft, ispden, mu) + phim * v1r_qbz(1, ifft, ispden, mu)
             end do

           end do
         end do
       end do
      end if

   end do ! iqst

   ABI_FREE(v1r_qibz)
 end do ! iq_ibz
 ABI_CHECK(iqst == nqbz, "iqst /= nqbz")

 db%v1scf_rpt = db%v1scf_rpt / nqbz
 call xmpi_sum(db%v1scf_rpt, comm, ierr)

 ! Build mpi_type for executing fourdp in sequential.
 call ngfft_seq(ngfft_qspace, [nq1, nq2, nq3])
 call initmpi_seq(mpi_enreg_seq)
 call init_distribfft_seq(mpi_enreg_seq%distribfft,'c',ngfft_qspace(2),ngfft_qspace(3),'all')
 call init_distribfft_seq(mpi_enreg_seq%distribfft,'f',ngfft_qspace(2),ngfft_qspace(3),'all')

 !ABI_STAT_MALLOC(all_v1qr, (nqbz, 2, nfft, db%nspden, db%natom3), ierr)
 !ABI_CHECK(ierr==0, "oom in all_v1qr")
 !all_v1qr = zero

 cnt = 0
 do mu=1,db%natom3
   do ispden=1,db%nspden
     do ifft=1,nfft
       !cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
       !call fourdp(1,all_v1qr(:,:,ifft,ispden,mu),db%v1scf_rpt(:,ifft,ispden,mu),+1,&
       ! mpi_enreg_seq,nqbz,ngfft_qspace,paral_kgb0,tim_fourdp0)
     end do
   end do
 end do

 !ABI_FREE(all_v1qr)
 call destroy_mpi_enreg(mpi_enreg_seq)

 ABI_FREE(emiqr)
 ABI_FREE(qibz)
 ABI_FREE(wtq)
 ABI_FREE(qbz)
 ABI_FREE(indqq)
 ABI_FREE(iperm)
 ABI_FREE(bz2ibz_sort)
 ABI_FREE(v1r_qbz)

end subroutine dvdb_ftinterp_setup
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_ftinterp_qpt
!! NAME
!!  dvdb_ftinterp_qpt
!!
!! FUNCTION
!!  Fourier interpolation of potentials for a given q-point
!!
!! INPUTS
!!  qpt(3)=q-point in reduced coordinates.
!!  nfft=Number of fft-points treated by this processors
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  comm=MPI communicator
!!
!! OUTPUT
!!  ov1r(2*nfft, nspden, 3*natom)=Interpolated DFPT potentials at the given q-point.
!!
!! PARENTS
!!      m_phgamma,mrgdv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_ftinterp_qpt(db, qpt, nfft, ngfft, ov1r, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_ftinterp_qpt'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,comm
 type(dvdb_t),intent(in) :: db
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(out) :: ov1r(2,nfft,db%nspden,db%natom3)

!Local variables-------------------------------
!scalars
 integer,parameter :: cplex2=2
 integer :: ir,ispden,ifft,mu,idir,ipert,timerev_q,nproc,my_rank,cnt,ierr
 real(dp) :: wr,wi
!arrays
 integer :: symq(4,2,db%cryst%nsym)
 real(dp),allocatable :: eiqr(:,:)

! *************************************************************************

 ABI_CHECK(allocated(db%v1scf_rpt), "v1scf_rpt is not allocated")

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 ! Examine the symmetries of the q wavevector
 call littlegroup_q(db%cryst%nsym,qpt,symq,db%cryst%symrec,db%cryst%symafm,timerev_q,prtvol=db%prtvol)

 ! Compute FT phases for this q-point.
 ABI_MALLOC(eiqr, (2,db%nrpt))
 call calc_eiqr(qpt,db%nrpt,db%rpt,eiqr)

 ov1r = zero
 do mu=1,db%natom3
   idir = mod(mu-1, 3) + 1; ipert = (mu - idir) / 3 + 1
   do ispden=1,db%nspden
     do ifft=1,nfft
       ! MPI-parallelism
       cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle

       do ir=1,db%nrpt
         wr = db%v1scf_rpt(1, ir,ifft,ispden,mu)
         wi = db%v1scf_rpt(2, ir,ifft,ispden,mu)
         ov1r(1,ifft,ispden,mu) = ov1r(1,ifft,ispden,mu) + wr*eiqr(1,ir) - wi * eiqr(2,ir)
         ov1r(2,ifft,ispden,mu) = ov1r(2,ifft,ispden,mu) + wr*eiqr(2,ir) + wi * eiqr(1,ir)
       end do
     end do

     ! Remove the phase.
     call times_eikr(-qpt, ngfft, nfft, 1, ov1r(:,:,ispden,mu))
   end do

   ! Be careful with gamma and cplex!
   if (db%symv1) then
     call v1phq_symmetrize(db%cryst,idir,ipert,symq,ngfft,cplex2,nfft,db%nspden,db%nsppol,db%mpi_enreg,ov1r(:,:,:,mu))
   end if
 end do ! mu

 call xmpi_sum(ov1r, comm, ierr)

 ABI_FREE(eiqr)

end subroutine dvdb_ftinterp_qpt
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_findq
!! NAME
!!  dvdb_findq
!!
!! FUNCTION
!!  Find the index of the q-point in db%qpts. Non zero umklapp vectors are not allowed.
!!  Returns -1 if not found.
!!
!! INPUTS
!!  qpt(3)=q-point in reduced coordinates.
!!  [qtol]=Optional tolerance for q-point comparison.
!!         For each reduced direction the absolute difference between the coordinates must be less that qtol
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function dvdb_findq(db, qpt, qtol) result(iqpt)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_findq'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),optional,intent(in) :: qtol
 type(dvdb_t),intent(in) :: db
!arrays
 real(dp),intent(in) :: qpt(3)

!Local variables-------------------------------
!scalars
 integer :: iq
 real(dp) :: my_qtol

! *************************************************************************

 my_qtol = 0.0001_dp; if (present(qtol)) my_qtol = qtol

 iqpt = -1
 do iq=1,db%nqpt
   if (all(abs(db%qpts(:, iq) - qpt) < my_qtol)) then
      iqpt = iq; exit
   end if
 end do

end function dvdb_findq
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_seek
!! NAME
!!  dvdb_seek
!!
!! FUNCTION
!!   Move the internal file pointer so that it points to the
!!   block with (idir, ipert, iqpt). Needed only if dvdb%iomode==IO_MODE_FORTRAN
!!
!! INPUTS
!!   idir,ipert,iqpt = (direction, perturbation, q-point) indices
!!
!! SIDE EFFECTS
!!   db<type(dvdb_t)> : modifies db%f90_fptr and the internal F90 file pointer.
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_seek(db, idir, ipert, iqpt)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_seek'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in)  :: idir, ipert, iqpt
 type(dvdb_t),intent(inout) :: db

!Local variables-------------------------------
 integer :: pos_wanted,ii,ispden,ierr
 real(dp),parameter :: fake_qpt(3)=zero
 !character(len=500) :: msg

! *************************************************************************

 if (db%iomode == IO_MODE_FORTRAN) then
   ! Numb code: rewind the file and read it from the beginning
   ! TODO: Here I need a routine to backspace the hdr!
   call dvdb_rewind(db)
   pos_wanted = db%pos_dpq(idir,ipert,iqpt)

   do ii=1,pos_wanted-1
     !write(std_out,*)"in seek with ii: ",ii,"pos_wanted: ",pos_wanted
     ierr = my_hdr_skip(db%fh, -1, -1, fake_qpt)
     ABI_CHECK(ierr==0, "Error in hdr_skip")
     ! Skip the records with v1.
     do ispden=1,db%nspden
       read(db%fh)
     end do
   end do
   db%current_fpos = pos_wanted

 else
   MSG_ERROR("should not be called when iomode /= IO_MODE_FORTRAN")
 end if

end subroutine dvdb_seek
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_rewind
!! NAME
!!  dvdb_rewind
!!
!! FUNCTION
!!   Rewind the file and move to the first header. Needed only if dvdb%iomode==IO_MODE_FORTRAN
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_rewind(db)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_rewind'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(dvdb_t),intent(inout) :: db

! *************************************************************************

 if (db%iomode == IO_MODE_FORTRAN) then
   rewind(db%fh)
   read(db%fh)  ! version
   read(db%fh)  ! numv1
   db%current_fpos = 1

 else
   MSG_ERROR("should not be called when iomode /= IO_MODE_FORTRAN")
 end if

end subroutine dvdb_rewind
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/my_hdr_skip
!! NAME
!!  my_hdr_skip
!!
!! FUNCTION
!!  Skip the header without rewinding the file. Return exit code.
!!
!! NOTES
!!  Because hdr_skip rewinds the file and I'm not gonna change that ugly code.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function my_hdr_skip(unit,idir,ipert,qpt) result(ierr)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'my_hdr_skip'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unit,idir,ipert
 real(dp),intent(in) :: qpt(3)

!Local variables-------------------------------
 integer :: fform
 type(hdr_type) :: tmp_hdr
!************************************************************************

 ierr = 0
 call hdr_fort_read(tmp_hdr, unit, fform)
 ABI_CHECK(fform /= 0, "hdr_fort_read returned fform == 0")

 if (idir /= -1 .and. ipert /= -1) then
   if (idir /= mod(tmp_hdr%pertcase-1, 3) + 1 .or. &
       ipert /= (tmp_hdr%pertcase - idir) / 3 + 1 .or. &
       any(abs(qpt - tmp_hdr%qptn) > tol14)) ierr = -1
 end if

 call hdr_free(tmp_hdr)
 if (fform == 0) ierr = 1

end function my_hdr_skip
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_list_perts
!! NAME
!!  dvdb_list_perts
!!
!! FUNCTION
!!  Given a q-point mesh, this routine checks if all the (phonon) perturbations
!!  are available taking into account symmetries.
!!
!! INPUTS
!!  ngqpt(3)=Q-mesh divisions. If all(ngqpt == -1), the list of q-points in the DVDB
!!    (i.e. db%qpts) is analyzed instead of the q-points generated from ngqpt.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      eph,m_dvdb,mrgdv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_list_perts(db, ngqpt)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_list_perts'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(dvdb_t),target,intent(in) :: db
!arrays
 integer,intent(in) :: ngqpt(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: iout0=0,iscf1=1
 integer :: tot_miss,tot_weird,miss_q,idir,ipert,nqpt_computed,iv1,psy,weird_q
 integer :: iq_ibz,nqibz,iq_file,chksymbreak,kptopt,nshiftk,ii,timerev_q
 real(dp) :: kptrlen
 character(len=500) :: msg,ptype,found
 type(crystal_t),pointer :: cryst
!arrays
 integer,parameter :: vacuum0(3)=[0, 0, 0]
 integer :: rfdir(3),kptrlatt(3,3)
 integer,allocatable :: pertsy(:,:),symq(:,:,:),rfpert(:)
 real(dp) :: qq(3),shiftk(3,210)
 real(dp),allocatable :: qibz(:,:),wtq(:)

! *************************************************************************

 cryst => db%cryst

 if (all(ngqpt == -1)) then
   ! Will test the q-points in db
   call alloc_copy(db%qpts, qibz)
   nqibz = db%nqpt
 else
   ! Will test the q-points in the IBZ associated to ngqpt hence
   ! I call `getkgrid` to get the q-points in the IBZ from ngqpt.
   chksymbreak = 0; kptopt = 1; shiftk = zero; nshiftk = 1; kptrlatt = 0
   do ii=1,3
     kptrlatt(ii, ii) = ngqpt(ii)
   end do

   ABI_MALLOC(qibz, (3,0))
   ABI_MALLOC(wtq, (0))

   call getkgrid(chksymbreak,iout0,iscf1,qibz,kptopt,kptrlatt,kptrlen,&
     cryst%nsym,0,nqpt_computed,nshiftk,cryst%nsym,cryst%rprimd,shiftk,cryst%symafm,cryst%symrel,vacuum0,wtq)

   ABI_FREE(qibz)
   ABI_FREE(wtq)

   nqibz = nqpt_computed
   ABI_MALLOC(qibz, (3,nqibz))
   ABI_MALLOC(wtq, (nqibz))

   call getkgrid(chksymbreak,iout0,iscf1,qibz,kptopt,kptrlatt,kptrlen,&
     cryst%nsym,nqibz,nqpt_computed,nshiftk,cryst%nsym,cryst%rprimd,shiftk,cryst%symafm,cryst%symrel,vacuum0,wtq)

   ABI_FREE(wtq)
 end if

 ! Initialize the list of perturbations rfpert and rdfir
 ! WARNING: Only phonon perturbations are considered for the time being.
 ABI_MALLOC(rfpert,(db%mpert))
 rfpert = 0; rfpert(1:cryst%natom) = 1; rfdir = 1

 ABI_MALLOC(symq, (4,2,cryst%nsym))
 ABI_MALLOC(pertsy, (3,db%mpert))

 ! Loop over the q-points in the IBZ and test whether the q-point is present
 ! and if all the independent perturbations are available.
 ! `tot_miss` counts the number of irreducible perturbations not found in the DVDB (critical)
 ! `tot_weird` counts the number of independent perturbations found in the DVDB (not critical)
 tot_miss = 0; tot_weird = 0
 do iq_ibz=1,nqibz
   qq = qibz(:,iq_ibz)
   iq_file = dvdb_findq(db, qq)

   ! Examine the symmetries of the q wavevector
   call littlegroup_q(cryst%nsym,qq,symq,cryst%symrec,cryst%symafm,timerev_q,prtvol=db%prtvol)

   ! Determine the symmetrical perturbations. Meaning of pertsy:
   !    0 for non-target perturbations
   !    1 for basis perturbations
   !   -1 for perturbations that can be found from basis perturbations
   call irreducible_set_pert(cryst%indsym,db%mpert,cryst%natom,cryst%nsym,&
     pertsy,rfdir,rfpert,symq,cryst%symrec,cryst%symrel)

   if (iq_file /= -1) then
     ! This q-point is in the db. Test if all the independent perturbations are available.
     msg = sjoin("qpoint:",ktoa(qq),"is present in the DVDB file")
     call wrtout(std_out,msg,"COLL")
     call wrtout(std_out,' The list of irreducible perturbations for this q vector is:','COLL')
     ii = 0; weird_q = 0; miss_q = 0
     do ipert=1,db%mpert
       do idir=1,3
         psy = pertsy(idir,ipert)
         if (psy == 0) cycle
         iv1 = db%pos_dpq(idir,ipert,iq_file)
         ptype = "independent"; if (psy == -1) ptype = "symmetric"
         found = "Yes"; if (iv1 == 0) found = "No"

         if (psy == 1 .and. iv1 == 0) miss_q = miss_q + 1
         if (psy == -1 .and. iv1 /= 0) weird_q = weird_q + 1

         ii=ii+1
         write(msg,'(i5,a,i2,a,i4,4a)')ii,')  idir=',idir,', ipert=',ipert,", type=",trim(ptype),", found=",trim(found)
         call wrtout(std_out,msg,'COLL')
       end do
     end do

     if (weird_q /= 0) then
       write(msg,"(a,i0,a)")"DVDB is overcomplete. ",weird_q, " perturbation(s) can be reconstructed by symmetry."
       call wrtout(std_out,msg,"COLL")
     end if

     tot_weird = tot_weird + weird_q
     tot_miss = tot_miss + miss_q
     if (miss_q /=0) then
       write(msg,"(a,i0,a)")"WARNING: ",miss_q, " independent perturbation(s) are missing!. "
       call wrtout(std_out,msg,"COLL")
     end if

   else
     ! This q-point is not present in dvdb. Print the list of independent perturbations.
     msg = sjoin("qpoint:",ktoa(qq),"is NOT present in the DVDB file")
     call wrtout(std_out,msg,"COLL")
     call wrtout(std_out,' The list of irreducible perturbations for this q vector is:','COLL')
     ii = 0
     do ipert=1,db%mpert
       do idir=1,3
         if (pertsy(idir,ipert) == 1) then
           ii=ii+1
           write(msg,'(i5,a,i2,a,i4,a)')ii,')  idir=',idir,', ipert=',ipert,", type=independent, found=No"
           call wrtout(std_out,msg,'COLL')
           tot_miss = tot_miss + 1
         end if
       end do
     end do
   end if

   call wrtout(std_out,"","COLL")
 end do ! iq_ibz

 if (tot_miss /= 0) then
    write(msg,"(2a,i0,a)")ch10,"There are ",tot_miss, " independent perturbations missing!. "
    call wrtout(std_out,msg,"COLL")
 else
    call wrtout(std_out,"All the independent perturbations are available","COLL")
    if (tot_weird /= 0) then
      call wrtout(std_out,"Note however that the DVDB is overcomplete as symmetric perturbations are present.", "COLL")
    end if
 end if

 ABI_FREE(qibz)
 ABI_FREE(rfpert)
 ABI_FREE(symq)
 ABI_FREE(pertsy)

end subroutine dvdb_list_perts
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_check_v1sym
!! NAME
!!  dvdb_check_v1sym
!!
!! FUNCTION
!!  Debugging tool used to check whether the DFPT potentials in real space fulfill
!!  the correct symmetries on the real space FFT mesh.
!!
!! INPUTS
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_check_v1sym(db)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_check_v1sym'
 use interfaces_32_util
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(dvdb_t),target,intent(inout) :: db

!Local variables-------------------------------
!scalars
 integer,parameter :: iout0=0,rfmeth2=2,syuse0=0
 integer :: iqpt,timerev_q,idir,ipert,nsym1,cplex,v1pos
 integer :: isym,nfft,ifft,ifft_rot,ispden
 real(dp) :: max_err,re,im,vre,vim !,pre,pim
 character(len=500) :: msg
 logical :: isok
 type(crystal_t),pointer :: cryst
!arrays
 integer :: ngfft(18)
 integer,allocatable :: symq(:,:,:),symafm1(:),symrel1(:,:,:),irottb(:,:)
 real(dp) :: qpt(3)
 real(dp),allocatable :: tnons1(:,:),v1scf(:,:)

! *************************************************************************

 ABI_CHECK(db%nspinor==1, "nspinor == 2 not coded")

 cryst => db%cryst
 ABI_MALLOC(symq, (4,2,cryst%nsym))
 ABI_MALLOC(symafm1, (cryst%nsym))
 ABI_MALLOC(symrel1, (3,3,cryst%nsym))
 ABI_MALLOC(tnons1, (3,cryst%nsym))

 do iqpt=1,db%nqpt
   qpt = db%qpts(:, iqpt)

   ! Examine the symmetries of the q wavevector
   call littlegroup_q(cryst%nsym,qpt,symq,cryst%symrec,cryst%symafm,timerev_q,prtvol=db%prtvol)

   do ipert=1,db%natom
     do idir=1,3
       v1pos = db%pos_dpq(idir, ipert, iqpt); if (v1pos == 0) cycle

       ! Determines the set of symmetries that leaves the perturbation invariant.
       call littlegroup_pert(cryst%gprimd,idir,cryst%indsym,std_out,ipert,cryst%natom,cryst%nsym,nsym1,rfmeth2,&
         cryst%symafm,symafm1,symq,cryst%symrec,cryst%symrel,symrel1,syuse0,cryst%tnons,tnons1)

       ! FIXME: (3) or (18)
       cplex = db%cplex_v1(v1pos)
       ngfft(1:3) = db%ngfft3_v1(:, v1pos)
       nfft = product(ngfft(:3))
       ABI_MALLOC(v1scf, (cplex*nfft, db%nspden))

       if (dvdb_read_onev1(db, idir, ipert, iqpt, cplex, nfft, ngfft, v1scf, msg) /= 0) then
         MSG_ERROR(msg)
       end if

       ABI_MALLOC(irottb, (nfft,nsym1))
       call rotate_fft_mesh(nsym1,symrel1,tnons1,ngfft,irottb,isok)
       if (.not. isok) then
         MSG_WARNING("Real space FFT mesh is not compatible with symmetries!")
       end if

       max_err = zero
       do isym=1,nsym1
         do ispden=1,db%nspden
           do ifft=1,nfft
             ifft_rot = irottb(ifft, isym)
             !pre =  cos(two_pi * dot_product(qpt, tnons1(:,isym)))
             !pim = -sin(two_pi * dot_product(qpt, tnons1(:,isym)))
             if (cplex == 2) then
               re = v1scf(2*ifft_rot-1, ispden)
               im = v1scf(2*ifft_rot  , ispden)
               !vre = re * pre - im * pim
               !vim = re * pim + im * pre
               vre = re; vim = im

               re = v1scf(2*ifft-1, ispden) - vre
               im = v1scf(2*ifft  , ispden) - vim
             else
               re = v1scf(ifft, ispden) - v1scf(ifft_rot, ispden)
               im = zero
             end if
             !if (sqrt(re**2 + im**2) > tol6) write(std_out,*)"ifft,isym,err: ",ifft,isym,sqrt(re**2 + im**2)
             max_err = max(max_err, sqrt(re**2 + im**2))
           end do
         end do
       end do
       write(std_out,"(3(a,i0),a,es16.8)")"For iqpt= ",iqpt,", idir= ",idir,", ipert= ",ipert,",max_err= ",max_err

       ABI_FREE(irottb)
       ABI_FREE(v1scf)
     end do
   end do

 end do ! iqpt

 ABI_FREE(symq)
 ABI_FREE(symafm1)
 ABI_FREE(symrel1)
 ABI_FREE(tnons1)

end subroutine dvdb_check_v1sym
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_test_symmetries
!! NAME
!!  dvdb_test_symmetries
!!
!! FUNCTION
!!  Debugging tool used to test the symmetrization of the DFPT potentials.
!!
!! INPUTS
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      mrgdv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_test_symmetries(db)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_test_symmetries'
 use interfaces_32_util
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(dvdb_t),target,intent(inout) :: db

!Local variables-------------------------------
!scalars
 integer :: iqpt,pcase,idir,ipert,cplex,nfft,ispden,timerev_q !,ifft
 !character(len=500) :: msg
 type(crystal_t),pointer :: cryst
!arrays
 integer :: ngfft(18),rfdir(3),pflag(3, db%natom)
 real(dp) :: qpt(3)
 integer,allocatable :: pertsy(:,:),rfpert(:),symq(:,:,:)
 real(dp),allocatable :: file_v1scf(:,:,:,:),symm_v1scf(:,:,:,:)

! *************************************************************************

 cryst => db%cryst
 call ngfft_seq(ngfft, db%ngfft3_v1(:,1))
 nfft = product(ngfft(:3))

 ! Initialize the list of perturbations rfpert and rdfir
 ! WARNING: Only phonon perturbations are considered for the time being.
 ABI_MALLOC(rfpert,(db%mpert))
 rfpert = 0; rfpert(1:cryst%natom) = 1; rfdir = 1

 ABI_MALLOC(symq, (4,2,cryst%nsym))
 ABI_MALLOC(pertsy, (3,db%mpert))

 do iqpt=1,db%nqpt
   qpt = db%qpts(:,iqpt)

   ! Examine the symmetries of the q wavevector
   call littlegroup_q(cryst%nsym,qpt,symq,cryst%symrec,cryst%symafm,timerev_q,prtvol=db%prtvol)

   ! Determine the symmetrical perturbations. Meaning of pertsy:
   !    0 for non-target perturbations
   !    1 for basis perturbations
   !   -1 for perturbations that can be found from basis perturbations
   call irreducible_set_pert(cryst%indsym,db%mpert,cryst%natom,cryst%nsym,&
     pertsy,rfdir,rfpert,symq,cryst%symrec,cryst%symrel)

   ! Read all potentials (here I assume that all perturbations are available)
   call dvdb_readsym_allv1(db, iqpt, cplex, nfft, ngfft, file_v1scf, db%comm)

   ! Copy basis perturbations in symm_v1scf and set pflag
   ABI_MALLOC(symm_v1scf, (cplex,nfft,db%nspden,db%natom3))
   symm_v1scf = huge(one); pflag = 0
   do pcase=1,3*db%cryst%natom
     idir = mod(pcase-1, 3) + 1; ipert = (pcase - idir) / 3 + 1
     if (pertsy(idir, ipert) == 1) then
       symm_v1scf(:,:,:,pcase) = file_v1scf(:,:,:,pcase)
       pflag(idir,ipert) = 1
     end if
   end do

   ! Complete potentials
   call v1phq_complete(cryst,qpt,ngfft,cplex,nfft,db%nspden,db%nsppol,db%mpi_enreg,db%symv1,pflag,symm_v1scf)

   ! Compare values.
   do pcase=1,3*db%cryst%natom
     idir = mod(pcase-1, 3) + 1; ipert = (pcase - idir) / 3 + 1
     if (pflag(idir,ipert) == 1) cycle

     write(std_out,"(2(a,i0),2a)")"For idir: ",idir, ", ipert: ", ipert, ", qpt: ",trim(ktoa(qpt))
     do ispden=1,db%nspden
       write(std_out,*)"  max(abs(f1-f2))",maxval(abs(file_v1scf(:,:,ispden,pcase) - symm_v1scf(:,:,ispden,pcase)))
       call vdiff_print(vdiff_eval(cplex,nfft,file_v1scf(:,:,ispden,pcase),symm_v1scf(:,:,ispden,pcase),cryst%ucvol))

       !if (cplex == 1) then
       !  do ifft=1,nfft
       !    write(std_out,*)file_v1scf(1,ifft,ispden,pcase),symm_v1scf(1,ifft,ispden,pcase)
       !  end do
       !else
       !  do ifft=1,nfft
       !    write(std_out,*)file_v1scf(1,ifft,ispden,pcase),symm_v1scf(1,ifft,ispden,pcase),&
       !                    file_v1scf(2,ifft,ispden,pcase),symm_v1scf(2,ifft,ispden,pcase)
       !  end do
       !end if
     end do
     write(std_out,*)""
   end do

   ABI_FREE(symm_v1scf)
   ABI_FREE(file_v1scf)
 end do

 ABI_FREE(rfpert)
 ABI_FREE(symq)
 ABI_FREE(pertsy)

end subroutine dvdb_test_symmetries
!!***

!----------------------------------------------------------------------

!!****f* m_dvdb/dvdb_merge_files
!! NAME
!!  dvdb_merge_files
!!
!! FUNCTION
!!  Merge a list of POT1 files.
!!
!! INPUT
!!  nfiles=Number of files to be merged.
!!  dvdb_path=Name of output DVDB file.
!!  prtvol=Verbosity level.
!!
!! SIDE EFFECTS
!!   v1files=List of file names to merge. This list could be changed if POT1 files in netcdf format is found.
!!
!! PARENTS
!!      mrgdv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvdb_merge_files(nfiles, v1files, dvdb_path, prtvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dvdb_merge_files'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfiles,prtvol
 character(len=*),intent(in) :: dvdb_path
 character(len=*),intent(inout) :: v1files(nfiles)

!Local variables-------------------------------
!scalars
 integer :: fform_pot=102
 integer :: ii,jj,fform,ount,cplex,nfft,ifft,ispden,nperts
 integer :: n1,n2,n3,v1_varid,ierr
 logical :: qeq0
 character(len=500) :: msg
 type(hdr_type),pointer :: hdr1
 type(dvdb_t) :: dvdb
!arrays
 integer :: units(nfiles)
 real(dp),allocatable :: v1(:)
 type(hdr_type),target,allocatable :: hdr1_list(:)

!************************************************************************

 if (file_exists(dvdb_path)) then
   MSG_ERROR(sjoin("Cannot overwrite existing file:", dvdb_path))
 end if

 ! If a file is not found, try the netcdf version and change v1files accordingly.
 do ii=1,nfiles
   if (nctk_try_fort_or_ncfile(v1files(ii), msg) /= 0) then
     MSG_ERROR(msg)
   end if
 end do

 ! Read the headers
 ABI_DT_MALLOC(hdr1_list, (nfiles))
 do ii=1,nfiles
   write(std_out,"(a,i0,2a)")"- Reading header of file [",ii,"]: ",trim(v1files(ii))

   if (endswith(v1files(ii), ".nc")) then
#ifdef HAVE_NETCDF
      NCF_CHECK(nctk_open_read(units(ii), v1files(ii), xmpi_comm_self))
      call hdr_ncread(hdr1_list(ii),units(ii),fform)
#endif
   else
     if (open_file(v1files(ii), msg, newunit=units(ii), form="unformatted", action="read", status="old") /= 0) then
       MSG_ERROR(msg)
     end if
     call hdr_fort_read(hdr1_list(ii), units(ii), fform)
     ABI_CHECK(fform /= 0, sjoin("fform=0 while reading:", v1files(ii)))
   end if
   if (prtvol > 0) call hdr_echo(hdr1_list(ii), fform, 3, unit=std_out)
   ! FIXME: fform is zero in POT files!
   !ABI_CHECK(fform /= 0, sjoin("fform=0 while reading: ", v1files(ii)))
   !write(std_out,*)"done", trim(v1files(ii))
   if (hdr1_list(ii)%pertcase == 0) then
     MSG_ERROR(sjoin("Found GS potential:", v1files(ii)))
   end if
 end do

 ! Validate headers.
 ! TODO: Should perform consistency check on the headers and rearrange them in blocks of q-points.
 nperts = size(hdr1_list)

 ! Write dvdb file (we only support fortran binary format)
 if (open_file(dvdb_path, msg, newunit=ount, form="unformatted", action="write", status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if
 write(ount)dvdb_last_version
 write(ount)nperts

 do ii=1,nperts
   jj = ii
   hdr1 => hdr1_list(jj)
   call hdr_fort_write(hdr1, ount, fform_pot, ierr)
   ABI_CHECK(ierr == 0, "hdr_fort_write returned ierr = 0")

   qeq0 = (hdr1%qptn(1)**2+hdr1%qptn(2)**2+hdr1%qptn(3)**2<1.d-14)
   cplex = 2; if (qeq0) cplex = 1
   nfft = product(hdr1%ngfft(1:3))
   n1 = hdr1%ngfft(1); n2 = hdr1%ngfft(2); n3 = hdr1%ngfft(3)

   ABI_MALLOC(v1, (cplex*nfft))

   if (.not. endswith(v1files(ii), ".nc")) then
      ! Fortran IO
      do ispden=1,hdr1%nspden
        read(units(jj)) (v1(ifft), ifft=1,cplex*nfft)
        write(ount) (v1(ifft), ifft=1,cplex*nfft)
      end do
   else
#ifdef HAVE_NETCDF
      ! Netcdf IO
      ! netcdf array has shape [cplex, n1, n2, n3, nspden]
      NCF_CHECK(nf90_inq_varid(units(ii), "first_order_potential", v1_varid))
      do ispden=1,hdr1%nspden
        NCF_CHECK(nf90_get_var(units(ii), v1_varid, v1, start=[1,1,1,ispden], count=[cplex, n1, n2, n3, 1]))
        write(ount) (v1(ifft), ifft=1,cplex*nfft)
      end do
#endif
   end if

   ABI_FREE(v1)
 end do

 close(ount)

 do ii=1,size(hdr1_list)
   call hdr_free(hdr1_list(ii))
   if (.not. endswith(v1files(ii), ".nc")) then
     close(units(ii))
   else
#ifdef HAVE_NETCDF
     NCF_CHECK(nf90_close(units(ii)))
#endif
   endif
 end do
 ABI_DT_FREE(hdr1_list)

 write(std_out,"(a,i0,a)")"Merged successfully ",nfiles," files"

 ! List available perturbations.
 call dvdb_init(dvdb, dvdb_path, xmpi_comm_self)
 call dvdb_print(dvdb)
 call dvdb_list_perts(dvdb, [-1,-1,-1])
 call dvdb_free(dvdb)

end subroutine dvdb_merge_files
!!***

!!****f* m_dvdb/calc_eiqr
!! NAME
!!  calc_eiqr
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

subroutine calc_eiqr(qpt,nrpt,rpt,eiqr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_eiqr'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nrpt
!arrays
 real(dp),intent(in) :: qpt(3),rpt(3,nrpt)
 real(dp),intent(out) :: eiqr(2,nrpt)

!Local variables -------------------------
!scalars
 integer :: ir
 real(dp) :: qr

! *********************************************************************

 do ir=1,nrpt
   qr = two_pi * dot_product(qpt, rpt(:,ir))
   eiqr(1,ir) = cos(qr); eiqr(2,ir) = sin(qr)
 end do

end subroutine calc_eiqr
!!***

END MODULE m_dvdb
