!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_nctk
!! NAME
!! m_nctk
!!
!! FUNCTION
!!  Tools and wrappers for NETCDF-IO.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!   Remove create_nc_file, write_var_netcdf, the output of OUT.nc is dangereous
!!     because we can create too many dimensions and get
!!    nf90_def_dim - NetCDF library returned:   NetCDF: NC_MAX_DIMS exceeded
!!   Moreover the multiple calls to redef render the IO very inefficient
!!   That part should be rationalized!
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_nctk

 use defs_basis
 use m_abicore
 use m_build_info
 use m_errors
 use iso_c_binding
 use m_xmpi
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_fstrings,  only : itoa, sjoin, lstrip, char_count, strcat, endswith, startswith, ltoa
 use m_io_tools,  only : pick_aname, delete_file, file_exists
 use m_yaml,      only : DTSET_IDX

 implicit none

 private
!!***

 integer,public,parameter :: nctk_noid = -huge(1)
 ! This value is used to signal to procedures that IO should not be performed.

 ! Basic variables
 character(len=*),public,parameter :: etsfio_file_format = "ETSF Nanoquanta"
 character(len=*),public,parameter :: etsfio_conventions = "http://www.etsf.eu/fileformats/"

 integer,public,parameter :: etsfio_charlen = 80
 ! The value corresponding to character_string_len

 real,public,parameter :: etsfio_version = 3.3
 ! This is clearly wrong because one should use strings instad of floats that cannot be represented exactly
 ! but, unfortunately, it's in the specifications and we have to live with it!

#ifdef HAVE_NETCDF
 integer,public,parameter :: nctk_max_dims = NF90_MAX_DIMS
 ! Maximum number of dimensions

 integer,public,parameter :: nctk_slen = NF90_MAX_NAME
 ! String length used for the names of dimensions and variables
 ! The maximum allowable number of characters

#else
 ! replacements
 integer,public,parameter :: nctk_max_dims = 7
 integer,public,parameter :: nctk_slen=256
#endif

 character(len=5),private,parameter :: NCTK_IMPLICIT_DIMS(10) = [ &
   "one  ", "two  ", "three", "four ", "five ", "six  ", "seven", "eight", "nine ", "ten  "]

!!****t* m_nctk/nctkerr_t
!! NAME
!! nctkerr_t
!!
!! FUNCTION
!!
!! SOURCE
!!
 type,private :: nctkerr_t
#ifdef HAVE_NETCDF
   integer :: ncerr = nf90_noerr
#else
   integer :: ncerr = 0
#endif
   integer :: line = 0
   character(len=fnlen) :: file = "Dummy File"
   character(len=2048) :: msg="No error detected"
 end type nctkerr_t
!!***

 type(nctkerr_t),private,save :: einfo

!!****t* m_nctk/ncfdim_t
!! NAME
!! nctkdim_t
!!
!! FUNCTION
!!  Stores the name and the value of a netcdf dimension
!!
!! SOURCE

 type,public :: nctkdim_t
   character(len=nctk_slen) :: name   ! name of the dimension.
   integer :: value                   ! value of the dimension.
   !integer :: id
 end type nctkdim_t
!!***

!!****t* m_nctk/nctkarr_t
!! NAME
!! nctkarr_t
!!
!! FUNCTION
!!  Stores the name and the shape of a netcdf array
!!
!! SOURCE

 type,public :: nctkarr_t
   character(len=nctk_slen) :: name        ! name of the array.
   character(len=4) :: dtype               ! string specifying the type.
   character(len=nctk_slen) :: shape_str   ! string with the shape. e.g. "dim1, dim2" for [dim1, dim2] array.
 end type nctkarr_t
!!***

!!****s* m_nctk/nctkvar_t
!! NAME
!!  nctkvar_t
!!
!! FUNCTION
!!  This structure stores variable information, such as
!!  name, NetCDF id, type, shape and dimensions. It contains the following elements:
!!
!! SOURCE

 type nctkvar_t

   integer :: id
   ! the id used by NetCDF to access this variable.

   integer :: xtype
   ! the type of the variable

   integer :: ndims
   ! the number of dimensions (0 for scalar variable).

   integer :: natts
   ! The number of attributes associated to the variable

   character(len=nctk_slen) :: name
   ! the variable name.

   character(len=nctk_slen) :: dimnames(nctk_max_dims)
   ! the name corresponding to each dimension
   ! Only if array variable, only (1:ndims) are relevent

   integer :: dimids(nctk_max_dims) = -1
   ! The id of the dimensions. only (1:ndims) are relevent

   integer :: dimlens(nctk_max_dims) = 0
   ! the size for each dimension if array variable, only (1:ndims) are relevent

   !character(len=nctk_slen) :: dimnames(nctk_max_dims)
   !character(len=nctk_slen), pointer :: ncattnames(:)
   ! * ncattnames: the name corresponding to all associated attributes

 end type nctkvar_t
 !!***

 public :: nctk_idname              ! Return the nc identifier from the name of the variable.
 public :: nctk_ncify               ! Append ".nc" to ipath if ipath does not end with ".nc"
 public :: nctk_string_from_occopt  ! Return human-readable string with the smearing scheme.
 public :: nctk_fort_or_ncfile      ! Test wheter a path exists (fortran or nc file) and
                                    ! select iomode depending on file extension.
 public :: nctk_try_fort_or_ncfile  ! Return fortran or netcdf filename depending on the existence of the file.
 public :: nctk_test_mpiio          ! Test at run-time whether the netcdf library supports parallel IO.

#ifdef HAVE_NETCDF
 public :: ab_define_var            ! Helper function used to declare a netcdf variable.

 ! Helper functions
 public :: nctk_open_read           ! Open a file in read-only mode.
 public :: nctk_open_create         ! Create a netcdf file for modifications.
 public :: nctk_open_modify         ! Open an already existent file for modifications.
 public :: nctk_add_etsf_header     ! Add the ETSF-IO header.
 public :: nctk_set_defmode         ! Set the file in define mode (metadata)
 public :: nctk_set_datamode        ! Set the file in datamode (IO)
 public :: nctk_set_collective      ! Set collective access for a netcdf variable

 public :: nctk_def_dims            ! Define dimensions in a netcdf file.
 interface nctk_def_dims
   module procedure nctk_def_one_dim
   module procedure nctk_def_dim_list
 end interface nctk_def_dims

 public :: nctk_set_atomic_units    ! Set the value of the attributes "units" and "scale_to_atomic_units".
 public :: nctk_def_basedims        ! Define the basic dimensions used in ETSF-IO files.
 public :: nctk_def_scalars_type    ! Declare a list of scalars of the given type.
 public :: nctk_def_iscalars        ! Declare a list of integer scalars.
 public :: nctk_def_dpscalars       ! Declare a list of double precision scalars.
 public :: nctk_write_iscalars      ! Write a list of integer scalars.
 public :: nctk_write_dpscalars     ! Write a list of double precision scalars.
 public :: nctk_defnwrite_ivars     ! Declare and write a list of integer scalars.
 public :: nctk_defnwrite_dpvars    ! Declare and write a list of double precisions scalars.
 public :: nctk_write_ibz           ! Write k-points in the IBZ and corresponding weights.

 public :: nctk_def_one_array
 public :: nctk_def_arrays          ! Define netcdf arrays.

 interface nctk_def_arrays
   module procedure nctk_def_one_array
   module procedure nctk_def_array_list
 end interface nctk_def_arrays

 public :: nctk_get_dim

 public :: nctk_write_datar
 public :: nctk_read_datar
 public :: nctk_defwrite_nonana_terms  ! Write phonon frequencies and displacements for q-->0
                                       ! in the presence of non-analytical behaviour.
 public :: nctk_defwrite_nonana_raman_terms   ! Write raman susceptiblities for q-->0
 public :: nctk_defwrite_raman_terms   ! Write raman susceptiblities and frequencies for q=0

#endif
 public :: create_nc_file              ! FIXME: Deprecated
 public :: write_var_netcdf            ! FIXME: Deprecated
 public :: write_eig                   ! FIXME: Deprecated

 !integer,save ABI_PROTECTED, public :: nctk_cache_size = 32000000
 ! If the cache_size is provided when opening a netCDF-4/HDF5 file, it will be used instead
 ! of the default (32000000) as the size, in bytes, of the HDF5 chunk cache.

 !integer,save ABI_PROTECTED, public :: nctk_cache_nelems = 1000
 ! If cache_nelems is provided when opening a netCDF-4/HDF5 file, it will be used instead
 ! of the default (1000) as the maximum number of elements in the HDF5 chunk cache.

 !real,save ABI_PROTECTED, public :: nctk_cache_premtion = 0.75
 ! If cache_preemption is provided when opening a netCDF-4/HDF5 file, it will be used
 ! instead of the default (0.75) as the preemption value for the HDF5 chunk cache.

 logical, save ABI_PROTECTED, public :: nctk_has_mpiio = .False.
 ! This flag is set to true if the netcdf library supports parallel IO.
 ! Cannot use CPP flags because nf90_open_par and other similar functions are always
 ! exported by netcdf. As a consequence we have to check at run-time if we can
 ! perform parallel IO and we use nctk_has_mpiio to select the IO algorithms.

CONTAINS

!!****f* m_nctk/nctk_idname
!! NAME
!!  nctk_idname
!!
!! FUNCTION
!!  Return the nc identifier from the name of the variable
!!
!! PARENTS
!!
!! SOURCE

integer function nctk_idname(ncid, varname) result(varid)

!Arguments ------------------------------------
 integer,intent(in) :: ncid
 character(len=*),intent(in) :: varname

!Local variables-------------------------------
!scalars
 integer :: ncerr
 character(len=1000) :: msg

! *********************************************************************

#ifdef HAVE_NETCDF
 ncerr = nf90_inq_varid(ncid, varname, varid)

 if (ncerr /= nf90_noerr) then
   write(msg,'(5a)')&
     "NetCDF library returned: ",trim(nf90_strerror(ncerr)),ch10,&
     "while trying to get the ncid of variable: ",trim(varname)
   MSG_ERROR(msg)
 end if
#else
 MSG_ERROR("Netcdf support is not activated")
 write(std_out,*)ncid,varname
#endif

end function nctk_idname
!!***

!----------------------------------------------------------------------

!!****f* m_nctk/nctk_ncify
!! NAME
!!  nctk_ncify
!!
!! FUNCTION
!!  Append ".nc" to ipath if ipath does not end with ".nc"
!!
!! SOURCE

function nctk_ncify(ipath) result(opath)

 character(len=fnlen),intent(in) :: ipath
 character(len=fnlen) :: opath

! *********************************************************************

 if (.not. endswith(ipath, ".nc")) then
   opath = trim(ipath)//'.nc'
 else
   opath = ipath
 end if

end function nctk_ncify
!!***


!----------------------------------------------------------------------

!!****f* m_nctk/nctk_string_from_occopt
!! NAME
!!  nctk_string_from_occopt
!!
!! FUNCTION
!!
!! SOURCE

pure function nctk_string_from_occopt(occopt) result(smearing)

 integer,intent(in) :: occopt
 character(len=etsfio_charlen) :: smearing

! *********************************************************************

 select case (occopt)
 case (3)
   smearing = "Fermi-Dirac"
 case (4)
   smearing = "cold smearing of N. Marzari with minimization of the bump"
 case (5)
   smearing = "cold smearing of N. Marzari with monotonic function in the tail"
 case (6)
   smearing = "Methfessel and Paxton"
 case (7)
   smearing = "gaussian"
 case (8)
   smearing = "uniform"
 case default
   smearing = "none"
 end select

end function nctk_string_from_occopt
!!***

!----------------------------------------------------------------------

!!****f* m_nctk/nctk_fort_or_ncfile
!! NAME
!!  nctk_fort_or_ncfile
!!
!! FUNCTION
!!  Return the iomode used to perform IO operations on filename.
!!  If filename does not exist, a similar file with extension `.nc` is tried
!!  and iomode is set to IO_MODE_ETSF if the file exists.
!!  This trick is used to run the Abinit test suite in netcdf mode without changing the input files.
!!  The modification (if any) is logged to std_out.
!!
!! SIDE EFFECTS
!!  filename=Tentative filename in input. Changed to netcdf file if input filename does not exist
!!   and a file with netcdf extension is found.
!!
!! OUTPUT
!!  iomode=Flag selecting the IO library. Set to IO_MODE_ETSF if netcdf file, else IO_MODE_MPI
!!    if MPI supports it, finally IO_MODE_FORTRAN
!!  errmsg=String with error message. Use `if (len_trim(errmsg) /= 0) MSG_ERROR(errmsg)`
!!    to handle possible errors in the caller.
!!
!! PARENTS
!!      conducti_nc,optic
!!
!! CHILDREN
!!
!! SOURCE

subroutine nctk_fort_or_ncfile(filename, iomode, errmsg)

 character(len=*),intent(inout) :: filename
 character(len=*),intent(out) :: errmsg
 integer,intent(out) :: iomode

! *********************************************************************
  errmsg = ""

 ! Default value
#ifdef HAVE_MPI_IO
 iomode = IO_MODE_MPI
#else
 iomode = IO_MODE_FORTRAN
#endif

 !  Checking the existence of data file
 if (.not.file_exists(filename)) then
   ! Trick needed to run Abinit test suite in netcdf mode.
   if (file_exists(nctk_ncify(filename))) then
     write(std_out,"(3a)")"- File: ",trim(filename)," does not exist but found netcdf file with similar name."
     filename = nctk_ncify(filename); iomode = IO_MODE_ETSF
   end if
   if (.not. file_exists(filename)) then
     errmsg = 'Missing file: '//trim(filename)
   end if
 end if

end subroutine nctk_fort_or_ncfile
!!***

!----------------------------------------------------------------------

!!****f* m_nctk/nctk_try_fort_or_ncfile
!! NAME
!!  nctk_try_fort_or_ncfile
!!
!! FUNCTION
!!  If filename does not exist, a similar file with extension `.nc` is tried
!!  This trick is used to run the Abinit test suite in netcdf mode without changing the input files.
!!  The modification (if any) is logged to unit (Default: std_out)
!!
!! SIDE EFFECTS
!!  filename=Tentative filename in input. Changed to netcdf file if input filename does not exist
!!   and a file with netcdf extension is found.
!!
!! OUTPUT
!!  errmsg=String with error message if return value /= 0
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function nctk_try_fort_or_ncfile(filename, errmsg, unit) result(ierr)

 character(len=*),intent(inout) :: filename
 character(len=*),intent(out) :: errmsg
 integer,optional,intent(in) :: unit

!Local variables-------------------------------
!scalars
 integer :: unt

! *********************************************************************

 unt = std_out; if (present(unit)) unt = unit
 ierr = 0; errmsg = ""

 if (.not.file_exists(filename)) then
   ! Try netcdf exists.
   if (file_exists(nctk_ncify(filename))) then
     if (unt /= dev_null) then
       write(unt,"(3a)")"- File: ",trim(filename)," does not exist but found netcdf file with similar name."
     end if
     filename = nctk_ncify(filename)
   end if
   if (.not. file_exists(filename)) then
     ierr = 1; errmsg = 'Missing file: '//trim(filename)
   end if
 end if

end function nctk_try_fort_or_ncfile
!!***

!----------------------------------------------------------------------

!!****f* m_nctk/nctk_test_mpiio
!! NAME
!!  nctk_test_mpiio
!!
!! FUNCTION
!!  Test at run-time whether the netcdf library supports parallel IO and
!!  set the value of the module variable `nctk_has_mpiio`.
!!  This is a COLLECTIVE routine that should be called by all processors
!!  in MPI_COMM_WORLD at the beginning of the calculation
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine nctk_test_mpiio()

!Local variables-------------------------------
!scalars
#ifdef HAVE_NETCDF_MPI
 integer,parameter :: master=0
 integer :: ierr,ncid,ncerr
 character(len=500) :: msg
 character(len=fnlen) :: apath
#endif

! *********************************************************************

 nctk_has_mpiio = .False.

#ifdef HAVE_NETCDF_MPI
 if (xmpi_comm_rank(xmpi_world) == master) then
   ! Try to open a file with hdf5.
   apath = pick_aname()
   ncerr = nf90_create(apath, cmode=ior(ior(nf90_netcdf4, nf90_mpiio), nf90_write), ncid=ncid, &
     comm=xmpi_comm_self, info=xmpio_info)

   if (ncerr == nf90_noerr) then
     nctk_has_mpiio = .True.
     call wrtout(std_out,"Netcdf library supports MPI-IO", "COLL")
   else if (ncerr == nf90_enopar) then
     ! This is the value returned by the C function ifndef USE_PARALLEL
     MSG_WARNING(sjoin("netcdf lib does not support MPI-IO and: ", nf90_strerror(ncerr)))
     nctk_has_mpiio = .False.
   else
     ! Maybe something wrong in the low-level layer!
     MSG_WARNING(sjoin("Strange, netcdf seems to support MPI-IO but: ", nf90_strerror(ncerr)))
     nctk_has_mpiio = .False.
   end if

   ncerr = nf90_close(ncid)
   call delete_file(apath, ierr)
 end if

 ! Master broadcast nctk_has_mpiio
 call xmpi_bcast(nctk_has_mpiio,master,xmpi_world,ierr)

 if (.not. nctk_has_mpiio) then
   write(msg,"(5a)")&
      "The netcdf library does not support parallel IO, see message above",ch10,&
      "Abinit won't be able to produce files in parallel e.g. when paral_kgb==1 is used.",ch10,&
      "Action: install a netcdf4+HDF5 library with MPI-IO support."
   MSG_WARNING(msg)
 end if
#endif

#ifdef HAVE_NETCDF_DEFAULT
 if (.not. nctk_has_mpiio) then
   MSG_ERROR("--netcdf-default is on but netcdf library does not support MPI-IO. Aborting now")
 end if
#endif

end subroutine nctk_test_mpiio
!!***

#ifdef HAVE_NETCDF

!!****f* m_nctk/str2xtype
!! NAME
!!  str2xtype
!!
!! FUNCTION
!!  Return the netcdf type from a string. Possible values:
!!    c or ch   for NF90_CHAR
!!    i or int  for NF90_INT
!!   sp         for NF90_FLOAT
!!   dp         for NF90_DOUBLE
!!
!! SOURCE

integer function str2xtype(string) result(xtype)

!Arguments ------------------------------------
 character(len=*),intent(in) :: string

! *********************************************************************

 !Type 	FORTRAN API Mnemonic 	Bits
 !byte      NF90_BYTE           8
 !char      NF90_CHAR           8
 !short     NF90_SHORT          16
 !int       NF90_INT            32
 !float     NF90_FLOAT          32
 !double    NF90_DOUBLE         64

 select case (string)
 case ("c", "ch", "char")
   xtype = NF90_CHAR
 case ("i", "int")
   xtype = NF90_INT
 case ("sp")
   xtype = NF90_FLOAT
 case ("dp")
   xtype = NF90_DOUBLE
 case default
   MSG_ERROR(sjoin("Invalid string type:", string))
 end select

end function str2xtype
!!***

!!****f* m_nctk/bail_if_ncerr
!! NAME
!!  bail_if_ncerr
!!
!! FUNCTION
!!
!! INPUTS
!!  ncerr=Netcdf error
!!  line=line number of the file where problem occurred
!!  file=name of the f90 file containing the caller
!!
!! SOURCE

logical function bail_if_ncerr(ncerr, file, line) result(bail)

!Arguments ------------------------------------
 integer,intent(in) :: ncerr
 character(len=*),optional,intent(in) :: file
 integer,optional,intent(in) :: line

! *********************************************************************

 bail = (ncerr /= nf90_noerr)

 if (bail) then
   einfo%ncerr = ncerr
   einfo%file = "Subroutine Unknown"; if (present(file)) einfo%file = file
   einfo%line = 0; if (present(line)) einfo%line = line
   ! Append Netcdf string to user-defined message.
   write(einfo%msg,'(2a)')'NetCDF library raised: ',trim(nf90_strerror(ncerr))
 end if

end function bail_if_ncerr
!!***

!----------------------------------------------------------------------

!!****f* m_nctk/nctk_open_read
!! NAME
!!  nctk_open_read
!!
!! FUNCTION
!!   Open a netcdf file in read-only mode. Return exit status.
!!
!! INPUTS
!!  ncid=Netcdf identifier.
!!  comm=MPI communicator.
!!
!! PARENTS
!!
!! SOURCE

integer function nctk_open_read(ncid, path, comm) result(ncerr)

!Arguments ------------------------------------
 integer,intent(out) :: ncid
 integer,intent(in) :: comm
 character(len=*),intent(in) :: path

!Local variables-------------------------------
 integer :: nprocs

! *********************************************************************
 nprocs = xmpi_comm_size(comm)

 ! Enforce netcdf4 only if the communicator contains more than one processor.
 if (nctk_has_mpiio .and. nprocs > 1) then
#ifdef HAVE_NETCDF_MPI
   ncerr = nf90_open(path, mode=ior(ior(nf90_netcdf4, nf90_mpiio), nf90_nowrite),&
                     comm=comm, info=xmpio_info, ncid=ncid)
#else
   ncerr = nf90_einval
   MSG_WARNING("Netcdf without MPI support. Cannot open file, will abort in caller")
#endif
   NCF_CHECK_MSG(ncerr, sjoin("opening file:", path))
 else
   ncerr = nf90_open(path, mode=nf90_nowrite, ncid=ncid)
   NCF_CHECK_MSG(ncerr, sjoin("opening file:", path))
   !if (ncerr /= NC_EHDFERR) then
   !  ncerr = nf90_open(path, mode=ior(ior(nf90_netcdf4), nf90_nowrite),&
   !                    comm=comm, info=xmpio_info, ncid=ncid)
   !end if
   if (nprocs > 1) then
     ncerr = nf90_einval
     MSG_WARNING("netcdf without MPI-IO support with nprocs > 1! Will abort in the caller")
   end if
 endif

end function nctk_open_read
!!***

!----------------------------------------------------------------------

!!****f* m_nctk/nctk_open_create
!! NAME
!!  nctk_open_create
!!
!! FUNCTION
!!  Create and open the netcdf file. Return exis status.
!!
!! INPUTS
!!  path=Name of the file
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  ncid=Netcdf identifier.
!!
!! PARENTS
!!
!! SOURCE

integer function nctk_open_create(ncid, path, comm) result(ncerr)

!Arguments ------------------------------------
 integer,intent(out) :: ncid
 integer,intent(in) :: comm
 character(len=*),intent(in) :: path

!Local variables-------------------------------
 integer :: input_len
 character(len=strlen) :: my_string

! *********************************************************************

 ! Always use mpiio mode (i.e. hdf5) if available so that one perform parallel parallel IO
 if (nctk_has_mpiio) then
   ncerr = nf90_einval
#ifdef HAVE_NETCDF_MPI
   call wrtout(std_out, sjoin("- Creating HDf5 file with MPI-IO support:", path))
   ! Believe it or not, I have to use xmpi_comm_self even in sequential to avoid weird SIGSEV in the MPI layer!
   ncerr = nf90_create(path, cmode=ior(ior(nf90_netcdf4, nf90_mpiio), nf90_write), ncid=ncid, &
     comm=comm, info=xmpio_info)
#endif
 else
   call wrtout(std_out, sjoin("- Creating netcdf file WITHOUT MPI-IO support:", path))
   ncerr = nf90_create(path, ior(nf90_clobber, nf90_write), ncid)
   if (xmpi_comm_size(comm) > 1) then
     MSG_WARNING("netcdf without MPI-IO support with nprocs > 1!")
   end if
 end if
 NCF_CHECK(ncerr)

 ! Write etsf_header: file format, version and conventions.
 NCF_CHECK(nf90_put_att(ncid, NF90_GLOBAL, "file_format", etsfio_file_format))
 NCF_CHECK(nf90_put_att(ncid, NF90_GLOBAL, "file_format_version", etsfio_version))
 NCF_CHECK(nf90_put_att(ncid, NF90_GLOBAL, "Conventions", etsfio_conventions))

 ! Add info on the code that produced this file. These are extensions not in the standard.
 NCF_CHECK(nf90_put_att(ncid, NF90_GLOBAL, "code", "Abinit"))
 NCF_CHECK(nf90_put_att(ncid, NF90_GLOBAL, "abinit_version", abinit_version))

 ! Define the basic dimensions used in ETSF-IO files.
 NCF_CHECK(nctk_def_basedims(ncid, defmode=.True.))

 if (len_trim(INPUT_STRING) /= 0) then
   ! Write string with input.
   my_string = trim(INPUT_STRING)
   if (DTSET_IDX /= -1 .and. index(INPUT_STRING, "jdtset ") == 0) then
     my_string = "jdtset " // itoa(DTSET_IDX) // ch10 // trim(INPUT_STRING)
   end if

   input_len = len_trim(my_string)
   NCF_CHECK(nctk_def_dims(ncid, nctkdim_t("input_len", input_len)))
   NCF_CHECK(nctk_def_arrays(ncid, nctkarr_t("input_string", "c", "input_len")))
   !print *, trim(INPUT_STRING)

   if (xmpi_comm_rank(comm) == 0) then
     NCF_CHECK(nctk_set_datamode(ncid))
     ! Pass my_string(1:input_len)) instead from trim(string) to avoid SIGSEV on higgs_intel_19.0_serial
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "input_string"), my_string(1:input_len)))
     NCF_CHECK(nctk_set_defmode(ncid))
   end if
 end if

end function nctk_open_create
!!***

!----------------------------------------------------------------------

!!****f* m_nctk/nctk_open_modify
!! NAME
!!  nctk_open_modfy
!!
!! FUNCTION
!!   Open an already existent netcdf file for modifications. Return exit status.
!!
!! INPUTS
!!  path=File name.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  ncid=Netcdf identifier.
!!
!! PARENTS
!!
!! SOURCE

integer function nctk_open_modify(ncid, path, comm) result(ncerr)

!Arguments ------------------------------------
 integer,intent(out) :: ncid
 integer,intent(in) :: comm
 character(len=*),intent(in) :: path

! *********************************************************************

 if (.not. nctk_has_mpiio .and. xmpi_comm_size(comm) > 1) then
   MSG_ERROR("netcdf without MPI-IO support with nprocs > 1!")
 end if

 if (xmpi_comm_size(comm) > 1 .or. nctk_has_mpiio) then
   call wrtout(std_out, sjoin("- Opening HDf5 file with MPI-IO support:", path))
#ifdef HAVE_NETCDF_MPI
   ncerr = nf90_open_par(path, cmode=ior(ior(nf90_netcdf4, nf90_mpiio), nf90_write), &
     comm=comm, info=xmpio_info, ncid=ncid)
   NCF_CHECK_MSG(ncerr, sjoin("nf90_open_par: ", path))
#else
   MSG_ERROR("nprocs > 1 but netcdf does not support MPI-IO")
#endif
 else
   call wrtout(std_out, sjoin("- Opening netcdf file without MPI-IO support:", path))
   ncerr = nf90_open(path, nf90_write, ncid)
   NCF_CHECK_MSG(ncerr, sjoin("nf90_open: ", path))
 end if

 ! Set file in define mode.
 NCF_CHECK(nctk_set_defmode(ncid))

end function nctk_open_modify
!!***

!----------------------------------------------------------------------

!!****f* m_nctk/nctk_add_etsf_header
!! NAME
!!  nctk_add_etsf_header
!!
!! FUNCTION
!!   Add the etsf-io header to a file associated to a netcdf file handler.
!!
!! INPUTS
!!  ncid=Netcdf file identifier.
!!  * version = the number of version to be created.
!!  * title = (optional) a title for the file (80 characters max).
!!  * history = (optional) the first line of history (1024 characters max).
!!  * with_etsf_header = (optional) if true, will create a header
!!                       as defined in the ETSF specifications (default is .true.).
!!                       When value is .false., arguments title, history and version
!!                       are ignored.
!!
!! PARENTS
!!
!! SOURCE

integer function nctk_add_etsf_header(ncid, title, history) result(ncerr)

!Arguments ------------------------------------
 integer,intent(in) :: ncid
 character(len=*),optional,intent(in) :: title,history

!Local variables-------------------------------
 !integer :: ncerr
 character(len=*), parameter :: file_format = "ETSF Nanoquanta"
 character(len=*), parameter :: conventions = "http://www.etsf.eu/fileformats/"
 real :: format_version = 3.3 ! Real is not a good choice for a version!

! *********************************************************************

 ncerr = nctk_set_defmode(ncid)
 if (ncerr /= nf90_noerr) return

 ! The file format
 ncerr = nf90_put_att(ncid, NF90_GLOBAL, "file_format", file_format)
 if (ncerr /= nf90_noerr) return

 ! The version
 ncerr = nf90_put_att(ncid, NF90_GLOBAL, "file_format_version", format_version)
 if (ncerr /= nf90_noerr) return

 ! The conventions
 ncerr = nf90_put_att(ncid, NF90_GLOBAL, "Conventions", conventions)
 if (ncerr /= nf90_noerr) return

 ! The history if present
 if (present(history)) then
   ncerr = nf90_put_att(ncid, NF90_GLOBAL, "history", history(1:min(1024, len(history))))
   if (ncerr /= nf90_noerr) return
 end if

 ! The title if present
 if (present(title)) then
   ncerr = nf90_put_att(ncid, NF90_GLOBAL, "title", title(1:min(80, len(title))))
   if (ncerr /= nf90_noerr) return
 end if

 ! These are extensions not in the standard.
 ! Add info on the code that produced this file
 ncerr = nf90_put_att(ncid, NF90_GLOBAL, "code", "Abinit")
 if (ncerr /= nf90_noerr) return

 ncerr = nf90_put_att(ncid, NF90_GLOBAL, "code_version", ABINIT_VERSION)
 if (ncerr /= nf90_noerr) return

end function nctk_add_etsf_header
!!***

!----------------------------------------------------------------------

!!****f* m_nctk/nctk_set_defmode
!! NAME
!!  nctk_set_defmode
!!
!! FUNCTION
!!   Set the file in define mode, return exit status.
!!
!! INPUTS
!!  ncid=Netcdf identifier.
!!
!! PARENTS
!!
!! SOURCE

integer function nctk_set_defmode(ncid) result(ncerr)

!Arguments ------------------------------------
 integer,intent(in) :: ncid

! *********************************************************************

 ncerr = nf90_redef(ncid)
 ! Use same trick as in etsf_io
 if (ncerr /= nf90_noerr .and. ncerr /= -39) then
   NCF_CHECK(ncerr)
 else
   ncerr = nf90_noerr
 end if

end function nctk_set_defmode
!!***

!----------------------------------------------------------------------

!!****f* m_nctk/nctk_set_datamode
!! NAME
!!  nctk_set_datamode
!!
!! FUNCTION
!!  Set the file in data mode. Return exit status
!!
!! INPUTS
!!  ncid=Netcdf identifier.
!!  [reserve]
!!
!! OUTPUT
!!  ncerr=Exit status
!!
!! PARENTS
!!
!! SOURCE

integer function nctk_set_datamode(ncid, reserve) result(ncerr)

!Arguments ------------------------------------
 integer,intent(in) :: ncid
 logical,intent(in),optional :: reserve

!Local variables-------------------------------
!scalars
 logical :: do_reserve

! *********************************************************************

 do_reserve = .False.; if (present(reserve)) do_reserve = reserve

 ncerr = nf90_enddef(ncid)

 ! Use same trick as in etsf_io
 ! neded otherwise netcdf complains if the file is already in def mode.
 if (ncerr /= nf90_noerr .and. ncerr /= -38) then
   NCF_CHECK(ncerr)
 else
   ncerr = nf90_noerr
 end if

 return

 ! TODO
 if (do_reserve) then
   ncerr = nf90_enddef(ncid)
   !ncerr = nf90__enddef(ncid)
 else
   ncerr = nf90_enddef(ncid)
 end if

end function nctk_set_datamode
!!***

!----------------------------------------------------------------------

!!****f* m_nctk/nctk_set_collective
!! NAME
!!  nctk_set_collective
!!
!! FUNCTION
!!  Use collective IO for the given netcdf variable. Return exit status.
!!
!! INPUTS
!!  ncid=Netcdf file identifier.
!!  varid=Netcdf variable identifier.
!!
!! PARENTS
!!
!! SOURCE

integer function nctk_set_collective(ncid, varid) result(ncerr)

!Arguments ------------------------------------
 integer,intent(in) :: ncid,varid

! *********************************************************************

  ncerr = nf90_einval
#ifdef HAVE_NETCDF_MPI
  ncerr = nf90_var_par_access(ncid, varid, nf90_collective)
#else
  MSG_ERROR("nctk_set_collective should not be called if NETCDF does not support MPI-IO")
  ABI_UNUSED((/ncid, varid/))
#endif

end function nctk_set_collective
!!***

!!****f* m_nctk/nctk_def_one_dim
!! NAME
!!  nctk_def_one_dim
!!
!! FUNCTION
!!  Define list of dimensions variables and write their values.
!!  Return immediately if error
!!
!! INPUTS
!!  ncid=Netcdf identifier.
!!  dimnames(:)=List of strings with the name of the dimensions
!!  values(:)=List of integer scalars
!!  [defmode]=If True, the nc file is set in define mode (default=False)
!!  [prefix]=Prefix added to varnames and dimensions. Empty string if not specified.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function nctk_def_one_dim(ncid, nctkdim, defmode, prefix) result(ncerr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 logical,optional,intent(in) :: defmode
 character(len=*),optional,intent(in) :: prefix
!arrays
 type(nctkdim_t),intent(in) :: nctkdim

!Local variables-------------------------------
 integer :: dimid,dimlen
 character(len=nctk_slen) :: dname
 character(len=500) :: msg
! *********************************************************************

 ncerr = nf90_noerr

 if (present(defmode)) then
   if (defmode) then
     NCF_CHECK(nctk_set_defmode(ncid))
   end if
 end if

 if (present(prefix)) then
   if (any(nctkdim%name == NCTK_IMPLICIT_DIMS)) then
     dname = nctkdim%name
   else
     dname = strcat(prefix, nctkdim%name)
   end if
 else
   dname = nctkdim%name
 end if

 ! if dimension already exists, test whether it has the same value else define new dim.
 ncerr = nf90_inq_dimid(ncid, dname, dimid)

 if (ncerr == nf90_noerr) then
   NCF_CHECK(nf90_inquire_dimension(ncid, dimid, len=dimlen))
   if (dimlen /= nctkdim%value) then
     write(msg, "(2a,2(a,i0))")&
        "dimension already exists with a different value",ch10,&
        "file = ", dimlen, "; write = ", nctkdim%value
     MSG_ERROR(msg)
   end if
 else
   ncerr = nf90_def_dim(ncid, dname, nctkdim%value, dimid)
   NCF_CHECK(ncerr)
 end if

end function nctk_def_one_dim
!!***

!----------------------------------------------------------------------

!!****f* m_nctk/nctk_def_dim_list
!! NAME
!!  nctk_def_dim_list
!!
!! FUNCTION
!!  Define list of dimensions variables and write their values.
!!  Return immediately if error
!!
!! INPUTS
!!  ncid=Netcdf identifier.
!!  dimnames(:)=List of strings with the name of the dimensions
!!  values(:)=List of integer scalars
!!  [defmode]=If True, the nc file is set in define mode (default=False)
!!  [prefix]=Prefix added to varnames and dimensions. Empty string if not specified.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function nctk_def_dim_list(ncid, nctkdims, defmode, prefix) result(ncerr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 logical,optional,intent(in) :: defmode
 character(len=*),optional,intent(in) :: prefix
!arrays
 type(nctkdim_t),intent(in) :: nctkdims(:)

!Local variables-------------------------------
!scalars
 integer :: ii

! *********************************************************************

 ncerr = nf90_noerr
 if (present(defmode)) then
   if (defmode) then
     NCF_CHECK(nctk_set_defmode(ncid))
   end if
 end if

 do ii=1,size(nctkdims)
   if (present(prefix)) then
     ncerr = nctk_def_one_dim(ncid, nctkdims(ii), prefix=prefix)
   else
     ncerr = nctk_def_one_dim(ncid, nctkdims(ii))
   end if
   if (ncerr /= nf90_noerr) return
 end do

end function nctk_def_dim_list
!!***

!----------------------------------------------------------------------

!!****f* m_nctk/nctk_set_atomic_units
!! NAME
!!  nctk_set_atomic_units
!!
!! FUNCTION
!!  Set the attributes "units" to "atomic units" and "scale_to_atomic_units" to one.
!!
!! INPUTS
!!  ncid=Netcdf identifier.
!!  varname=Name of the variable
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function nctk_set_atomic_units(ncid, varname) result(ncerr)

!Arguments ------------------------------------
 integer,intent(in) :: ncid
 character(len=*),intent(in) :: varname

!Local variables-------------------------------
!scalars
 integer :: varid

! *********************************************************************

 ncerr = nf90_noerr

 varid = nctk_idname(ncid, varname)
 NCF_CHECK(nf90_put_att(ncid, varid, "units", "atomic units"))
 NCF_CHECK(nf90_put_att(ncid, varid, "scale_to_atomic_units", one))

end function nctk_set_atomic_units
!!***

!----------------------------------------------------------------------

!!****f* m_nctk/nctk_def_basedims
!! NAME
!!  nctk_def_basedims
!!
!! FUNCTION
!!  Define the basic dimensions used in ETSF-IO files.
!!
!! INPUTS
!!  ncid=Netcdf identifier.
!!  [defmode]=If True, the nc file is set in define mode (default=False)
!!
!! PARENTS
!!      m_dfpt_io,m_dfptdb,m_header,m_phonons
!!
!! CHILDREN
!!
!! SOURCE

integer function nctk_def_basedims(ncid, defmode) result(ncerr)

!Arguments ------------------------------------
 integer,intent(in) :: ncid
 logical,optional,intent(in) :: defmode

! *********************************************************************

 ncerr = nf90_noerr

 if (present(defmode)) then
   if (defmode) then
     NCF_CHECK(nctk_set_defmode(ncid))
   end if
 end if

 ! Basic ETSF-IO dimensions that should be always present in the file.
 ncerr = nctk_def_dims(ncid, [&
   nctkdim_t("complex", 2), nctkdim_t("symbol_length", 2), nctkdim_t("character_string_length", etsfio_charlen),&
   nctkdim_t("number_of_cartesian_directions", 3), nctkdim_t("number_of_reduced_dimensions", 3),&
   nctkdim_t("number_of_vectors", 3)])
 NCF_CHECK(ncerr)

 ! Useful integers.
 ncerr = nctk_def_dims(ncid, [&
   nctkdim_t("one", 1), nctkdim_t("two", 2), nctkdim_t("three", 3), &
   nctkdim_t("four", 4), nctkdim_t("five", 5), nctkdim_t("six", 6), &
   nctkdim_t("seven", 7), nctkdim_t("eight", 8), nctkdim_t("nine", 9), nctkdim_t("ten", 10)])
 NCF_CHECK(ncerr)

end function nctk_def_basedims
!!***

!!****f* m_nctk/ab_define_var
!!
!! NAME
!! ab_define_var
!!
!! FUNCTION
!! Write the definition of a variable, including units and mnemonics
!!
!! INPUTS
!! ncid = Identifier of the netcdf dataset
!! var_dim_id = Identifier of the Dimensions
!! var_id     = Identifier of the variable
!! var_mnemo  = String of mnemonics
!! var_name   = String with the name of the variable
!! var_type   = NetCDF type of variable (NF90_DOUBLE, etc)
!! var_units  = String of units
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      m_abihist,m_bse_io,m_effective_potential,write_eig
!!
!! CHILDREN
!!
!! SOURCE

subroutine ab_define_var(ncid, var_dim_id, var_id, var_type, var_name, var_mnemo, var_units)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ncid
 integer, intent(out) :: var_id
 character(len=*), intent(in) :: var_mnemo,var_units,var_name
 integer,intent(in) :: var_type
!arrays
 integer,intent(in) :: var_dim_id(:)

!Local variables-------------------------------
!scalars
 integer :: ncerr

! *************************************************************************

#ifdef HAVE_NETCDF
 ncerr = nf90_def_var(ncid, trim(var_name), var_type, var_dim_id, var_id)
 NCF_CHECK_MSG(ncerr," define variable "//trim(var_name))

 ncerr = nf90_put_att(ncid, var_id,  "units",trim(var_units))
 NCF_CHECK_MSG(ncerr," define attribute for "//trim(var_name))

 ncerr = nf90_put_att(ncid, var_id,  "mnemonics", trim(var_mnemo))
 NCF_CHECK_MSG(ncerr," define attribute for "//trim(var_name))
#endif

end subroutine ab_define_var
!!***

!!****f* m_nctk/nctk_def_scalars_type
!! NAME
!!  nctk_def_scalars_type
!!
!! FUNCTION
!!  Define list of **scalar** variables with a given type, Return immediately if error
!!
!! INPUTS
!!  ncid=Netcdf identifier.
!!  varnames(:)=List of strings with the name of the variables
!!  xtype=Type of the variables
!!  [defmode]=If True, the nc file is set in define mode (default=False)
!!  [prefix]=Prefix added to varnames and dimensions. Empty string if not specified.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function nctk_def_scalars_type(ncid, varnames, xtype, defmode, prefix) result(ncerr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid,xtype
 logical,optional,intent(in) :: defmode
 character(len=*),optional,intent(in) :: prefix
!arrays
 character(len=*),intent(in) :: varnames(:)

!Local variables-------------------------------
!scalars
 integer :: ii,varid
 character(len=nctk_slen) :: vname
 type(nctkvar_t) :: var

! *********************************************************************

 ncerr = nf90_noerr
 if (present(defmode)) then
   if (defmode) then
     NCF_CHECK(nctk_set_defmode(ncid))
   end if
 end if

 ! Special case where dimension is null
 do ii=1,size(varnames)
   vname = varnames(ii)
   if (present(prefix)) vname = strcat(prefix, vname)

   if (nf90_inq_varid(ncid, vname, varid) == nf90_noerr)  then
       ! Variable already exists. Check if type and dimensions agree
       call var_from_id(ncid, varid, var)

       if (.not. (var%xtype == xtype .and. var%ndims == 0)) then
         MSG_ERROR("variable already exists with a different definition.")
       else
         cycle ! Dimension matches, skip definition.
       end if

   else
     ! Define variable since it doesn't exist.
     ncerr = nf90_def_var(ncid, vname, xtype, varid)
     NCF_CHECK(ncerr)
   end if
 end do

end function nctk_def_scalars_type
!!***

!!****f* m_nctk/nctk_def_iscalars
!! NAME
!!  nctk_def_iscalars
!!
!! FUNCTION
!!  Define list of integer **scalar** variables. Return immediately if error
!!
!! INPUTS
!!  ncid=Netcdf identifier.
!!  varnames(:)=List of strings with the name of the variables
!!  [defmode]=If True, the nc file is set in define mode (default=False)
!!  [prefix]=Prefix added to varnames and dimensions. Empty string if not specified.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function nctk_def_iscalars(ncid, varnames, defmode, prefix) result(ncerr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 logical,optional,intent(in) :: defmode
 character(len=*),optional,intent(in) :: prefix
!arrays
 character(len=*),intent(in) :: varnames(:)

! *********************************************************************

 if (present(defmode)) then
   ncerr = nctk_def_scalars_type(ncid, varnames, nf90_int, defmode=defmode)
 else
   if (present(prefix)) then
     ncerr = nctk_def_scalars_type(ncid, varnames, nf90_int, prefix=prefix)
   else
     ncerr = nctk_def_scalars_type(ncid, varnames, nf90_int)
   end if
 end if

end function nctk_def_iscalars
!!***

!!****f* m_nctk/nctk_def_dpscalars
!! NAME
!!  nctk_def_dpscalars
!!
!! FUNCTION
!!  Define list of double precision **scalar** variables. Return immediately if error
!!
!! INPUTS
!!  ncid=Netcdf identifier.
!!  varnames(:)=List of strings with the name of the variables
!!  [defmode]=If True, the nc file is set in define mode (default=False)
!!  [prefix]=Prefix added to varnames and dimensions. Empty string if not specified.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function nctk_def_dpscalars(ncid, varnames, defmode, prefix) result(ncerr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 logical,optional,intent(in) :: defmode
 character(len=*),optional,intent(in) :: prefix
!arrays
 character(len=*),intent(in) :: varnames(:)

!Local variables-------------------------------
 character(len=nctk_slen) :: prefix_

! *********************************************************************
 prefix_ = ""; if (present(prefix)) prefix_ = prefix

 if (present(defmode)) then
   ncerr = nctk_def_scalars_type(ncid, varnames, nf90_double, defmode=defmode, prefix=prefix_)
 else
   ncerr = nctk_def_scalars_type(ncid, varnames, nf90_double, prefix=prefix_)
 end if

end function nctk_def_dpscalars
!!***

!!****f* m_nctk/nctk_def_one_array
!! NAME
!!  nctk_def_one_array
!!
!! FUNCTION
!!  Define list of arrays with a given type, Return immediately if error
!!
!! INPUTS
!!  ncid=Netcdf identifier.
!!  nctk_array=Array descriptor.
!!  [defmode]=If True, the nc file is set in define mode (default=False)
!!  [prefix]=Prefix added to varnames and dimensions. Empty string if not specified.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE


integer function nctk_def_one_array(ncid, nctk_array, defmode, varid, prefix) result(ncerr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 logical,optional,intent(in) :: defmode
 integer,optional,intent(out) :: varid
 type(nctkarr_t),intent(in) :: nctk_array
 character(len=*),optional,intent(in) :: prefix

!Local variables-------------------------------
!scalars
 integer :: ii,xtype,prev,cnt,nn,vid
 character(len=500) :: msg
!arrays
 integer :: dimids(NF90_MAX_DIMS),dimvals(NF90_MAX_DIMS)
 character(len=nctk_slen) :: sarr(NF90_MAX_DIMS), string, pre, vname, dimname
 type(nctkvar_t) :: var

! *********************************************************************

 pre = ""; if (present(prefix)) pre = prefix

 ncerr = nf90_noerr
 if (present(defmode)) then
   if (defmode) then
     NCF_CHECK(nctk_set_defmode(ncid))
   end if
 end if

 xtype = str2xtype(nctk_array%dtype)
 vname = strcat(pre, nctk_array%name)

 ! Build array of strings with the dimensions.
 string = nctk_array%shape_str
 nn = char_count(string, ",")
 ABI_CHECK(nn <= NF90_MAX_DIMS, "Too many dimensions!")

 ! Parse dimension names and add prefix (if any).
 if (nn == 0) then
   cnt = 1
   dimname = lstrip(string)
   if (any(dimname == NCTK_IMPLICIT_DIMS)) then
     sarr(1) = dimname
   else
     sarr(1) = strcat(pre, dimname)
   end if
 else
   prev = 0; cnt = 0
   do ii=1,len_trim(string)
     if (string(ii:ii) == ",") then
       cnt = cnt + 1
       dimname = lstrip(string(prev+1:ii-1))
       if (any(dimname == NCTK_IMPLICIT_DIMS)) then
         sarr(cnt) = dimname
       else
         sarr(cnt) = strcat(pre, dimname)
       end if
       prev = ii
     end if
   end do
   cnt = cnt + 1
   dimname = lstrip(string(prev+1:ii-1))
   if (any(dimname == NCTK_IMPLICIT_DIMS)) then
     sarr(cnt) = dimname
   else
     sarr(cnt) = strcat(pre, dimname)
   end if
 end if

 ! Get dimids
 do ii=1,cnt
   NCF_CHECK_MSG(nf90_inq_dimid(ncid, sarr(ii), dimids(ii)), sarr(ii))
   NCF_CHECK(nf90_inquire_dimension(ncid, dimids(ii), len=dimvals(ii)))
 end do

 ! Check if dimension already exists.
 ! Variable already exists. Check if type and dimensions agree
 if (nf90_inq_varid(ncid, vname, vid) == nf90_noerr)  then
   call var_from_id(ncid, vid, var)
   if (.not. (var%xtype == xtype .and. var%ndims == cnt)) then
      write(msg,"(4a,2(2(a,i0),a))")&
        "variable ",trim(vname)," already exists with a different definition:",ch10,&
        "In file:     xtype = ",var%xtype,", ndims = ",var%ndims,ch10,&
        "From caller: xtype = ",xtype,", ndims = ",cnt,ch10
      MSG_ERROR(msg)
   end if
   if (any(dimvals(1:cnt) /= var%dimlens(1:var%ndims))) then
      write(msg,"(4a,2(3a))")&
        "variable ",trim(vname)," already exists but with different shape.",ch10,&
        "In file:     dims = ",trim(ltoa(var%dimlens(:var%ndims))),ch10,&
        "From caller  dims = ",trim(ltoa(dimvals(:cnt))),ch10
      MSG_ERROR(msg)
   end if
   if (present(varid)) varid = vid
   return
 end if

 ! Define variable since it doesn't exist.
 ncerr = nf90_def_var(ncid, vname, xtype, dimids(1:cnt), vid)
 NCF_CHECK(ncerr)

 if (present(varid)) varid = vid

end function nctk_def_one_array
!!***

!!****f* m_nctk/nctk_def_array_list
!! NAME
!!  nctk_def_array_list
!!
!! FUNCTION
!!  Define list of arrays with a given type, Return immediately if error
!!
!! INPUTS
!!  ncid=Netcdf identifier.
!!  nctk_arrays(:)=List of array descriptors.
!!  [defmode]=If True, the nc file is set in define mode (default=False)
!!  [prefix]=Prefix added to varnames and dimensions. Empty string if not specified.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE


integer function nctk_def_array_list(ncid, nctk_arrays, defmode, prefix) result(ncerr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 logical,optional,intent(in) :: defmode
 character(len=*),optional,intent(in) :: prefix
!arrays
 type(nctkarr_t),intent(in) :: nctk_arrays(:)

!Local variables-------------------------------
!scalars
 integer :: ia

! *********************************************************************

 ncerr = nf90_noerr
 if (present(defmode)) then
   if (defmode) then
     NCF_CHECK(nctk_set_defmode(ncid))
   end if
 end if

 do ia=1,size(nctk_arrays)
   if (present(prefix)) then
     NCF_CHECK(nctk_def_one_array(ncid, nctk_arrays(ia), prefix=prefix))
   else
     NCF_CHECK(nctk_def_one_array(ncid, nctk_arrays(ia)))
   end if
 end do

end function nctk_def_array_list
!!***

!!****f* m_nctk/nctk_write_iscalars
!! NAME
!!  nctk_write_iscalars
!!
!! FUNCTION
!!  Write a list of **scalar** integer variables. Return immediately if error
!!
!! INPUTS
!!  ncid=Netcdf identifier.
!!  varnames(:)=List of strings with the name of the variables
!!  values(:)=List of integer scalars
!!  [datamode]=If True, the nc file is set in data mode (default=False)
!!
!! OUTPUT
!!  ncerr=Exit status
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function nctk_write_iscalars(ncid, varnames, values, datamode) result(ncerr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 logical,optional,intent(in) :: datamode
!arrays
 integer,intent(in) :: values(:)
 character(len=*),intent(in) :: varnames(:)

!Local variables-------------------------------
!scalars
 integer :: ii,varid

! *********************************************************************

 ABI_CHECK(size(varnames) == size(values), "Different size in varnames, values")

 ncerr = nf90_noerr
 if (present(datamode)) then
   if (datamode) then
     NCF_CHECK(nctk_set_datamode(ncid))
   end if
 end if

 do ii=1,size(varnames)
   NCF_CHECK_MSG(nf90_inq_varid(ncid, varnames(ii), varid), sjoin("Inquiring: ", varnames(ii)))
   NCF_CHECK(nf90_put_var(ncid, varid, values(ii)))
 end do

end function nctk_write_iscalars
!!***

!!****f* m_nctk/nctk_write_dpscalars
!! NAME
!!  nctk_write_dpscalars
!!
!! FUNCTION
!!  Write a list of **scalar** real(dp) variables. Return immediately if error
!!
!! INPUTS
!!  ncid=Netcdf identifier.
!!  varnames(:)=List of strings with the name of the variables
!!  values(:)=List of real(dp) scalars
!!  [datamode]=If True, the nc file is set in data mode (default=False)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function nctk_write_dpscalars(ncid, varnames, values, datamode) result(ncerr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 logical,optional,intent(in) :: datamode
!arrays
 real(dp),intent(in) :: values(:)
 character(len=*),intent(in) :: varnames(:)

!Local variables-------------------------------
!scalars
 integer :: ii,varid

! *********************************************************************

 ncerr = nf90_noerr

 ABI_CHECK(size(varnames) == size(values), "Different size in varnames, values")

 if (present(datamode)) then
   if (datamode) then
     NCF_CHECK(nctk_set_datamode(ncid))
   end if
 end if

 do ii=1,size(varnames)
   NCF_CHECK(nf90_inq_varid(ncid, varnames(ii), varid))
   NCF_CHECK(nf90_put_var(ncid, varid, values(ii)))
 end do

end function nctk_write_dpscalars
!!***

!!****f* m_nctk/nctk_defnwrite_ivars
!! NAME
!!  nctk_defnwrite_ivars
!!
!! FUNCTION
!!  Define list of integer **scalar** variables and write their values.
!!  Return immediately if error
!!
!! INPUTS
!!  ncid=Netcdf identifier.
!!  varnames(:)=List of strings with the name of the variables
!!  values(:)=List of integer scalars
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function nctk_defnwrite_ivars(ncid, varnames, values) result(ncerr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
!arrays
 integer,intent(in) :: values(:)
 character(len=*),intent(in) :: varnames(:)

!Local variables-------------------------------
!scalars
 integer :: ii,varid

! *********************************************************************

 ABI_CHECK(size(varnames) == size(values), "Different size in varnames, values")

 ncerr = nctk_def_iscalars(ncid, varnames, defmode=.True.)
 NCF_CHECK(ncerr)

 NCF_CHECK(nctk_set_datamode(ncid))
 do ii=1,size(varnames)
   varid = nctk_idname(ncid, varnames(ii))
   NCF_CHECK(nf90_put_var(ncid, varid, values(ii)))
 end do

end function nctk_defnwrite_ivars
!!***

!!****f* m_nctk/nctk_defnwrite_dpvars
!! NAME
!!  nctk_defnwrite_dpvars
!!
!! FUNCTION
!!  Define list of real(dp) **scalar** variables and write their values.
!!  Return immediately if error
!!
!! INPUTS
!!  ncid=Netcdf identifier.
!!  varnames(:)=List of strings with the name of the variables
!!  values(:)=List of integer scalars
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function nctk_defnwrite_dpvars(ncid, varnames, values) result(ncerr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
!arrays
 real(dp),intent(in) :: values(:)
 character(len=*),intent(in) :: varnames(:)

!Local variables-------------------------------
!scalars
 integer :: ii,varid
!arrays

! *********************************************************************
 ncerr = nf90_noerr

 ABI_CHECK(size(varnames) == size(values), "Different size in varnames, values")

 ncerr = nctk_def_dpscalars(ncid, varnames, defmode=.True.)
 NCF_CHECK(ncerr)

 NCF_CHECK(nctk_set_datamode(ncid))
 do ii=1,size(varnames)
   !write(std_out,*)varnames(ii)
   varid = nctk_idname(ncid, varnames(ii))
   NCF_CHECK(nf90_put_var(ncid, varid, values(ii)))
 end do

end function nctk_defnwrite_dpvars
!!***

!!****f* m_nctk/nctk_write_ibz
!! NAME
!!  nctk_write_ibz
!!
!! FUNCTION
!!  Write the list of the k-points in the IBZ with the corresponding weights in Netcdf format.
!!  Mainly used for passing data to AbiPy. This routine should be called by master only.
!!
!! INPUTS
!!  fname=File name
!!  kpoints(:,:)=List of k-points
!!  weights(:)=K-point weights
!!
!! OUTPUT
!!  ncerr=Exit status
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function nctk_write_ibz(fname, kpoints, weights) result(ncerr)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: fname
!arrays
 real(dp),intent(in) :: kpoints(:,:),weights(:)

!Local variables-------------------------------
!scalars
 integer :: nkpts,ncid

! *********************************************************************

 ABI_CHECK(size(kpoints, dim=2) == size(weights), "size(kpoints, dim=2) != size(weights)")
 nkpts = size(kpoints, dim=2)

 NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), sjoin("Creating:", fname))

 ncerr = nctk_def_dims(ncid, [ &
   nctkdim_t("number_of_reduced_dimensions",3), nctkdim_t("number_of_kpoints", nkpts)], defmode=.True.)
 NCF_CHECK(ncerr)

 ncerr = nctk_def_array_list(ncid, [&
   nctkarr_t('reduced_coordinates_of_kpoints', "dp", "number_of_reduced_dimensions, number_of_kpoints"),&
   nctkarr_t('kpoint_weights', "dp", "number_of_kpoints")])
 NCF_CHECK(ncerr)

 NCF_CHECK(nctk_set_datamode(ncid))

 ncerr = nf90_put_var(ncid, nctk_idname(ncid, 'reduced_coordinates_of_kpoints'), kpoints)
 NCF_CHECK(ncerr)

 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'kpoint_weights'), weights))

 NCF_CHECK(nf90_close(ncid))

end function nctk_write_ibz
!!***

!!****f* m_nctk/nctk_get_dim
!! NAME
!!  nctk_get_dim
!!
!! FUNCTION
!!  Get the value of a dimension from its name.
!!
!! INPUTS
!!  ncid=Netcdf identifier.
!!  dimname=Name of the dimension.
!!  [datamode]=If True, the nc file is set in data mode (default=False)
!!
!! OUTPUT
!!  dimlen=Value of the dimension.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function nctk_get_dim(ncid, dimname, dimlen, datamode) result(ncerr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 character(len=*),intent(in) :: dimname
 integer,intent(out) :: dimlen
 logical,optional,intent(in) :: datamode

!Local variables-------------------------------
!scalars
 integer :: dimid

! *********************************************************************

 ncerr = nf90_noerr

 if (present(datamode)) then
   if (datamode) then
     NCF_CHECK(nctk_set_datamode(ncid))
   end if
 end if

 NCF_CHECK(nf90_inq_dimid(ncid, dimname, dimid))
 NCF_CHECK(nf90_inquire_dimension(ncid, dimid, len=dimlen))

end function nctk_get_dim
!!***

!----------------------------------------------------------------------

!!****f* m_nctk/nctk_write_datar
!! NAME
!! nctk_write_datar
!!
!! FUNCTION
!!  Write an array in real space in netcdf format
!!
!! INPUTS
!!  path=Filename
!!  varname=Name of the variable to write.
!!  ngfft(18)=information about 3D FFT
!!  cplex=1 for real arrays (e.g. GS rhor), 2 for complex array.
!!  nfft=number of points in the real space FFT mesh treated by this MPI proc
!!  nspden=number of spin-density components
!!  comm_fft=MPI communicator (used only if MPI-FFT).
!!  fftn3_distrib(n3)=rank of the processors which own fft planes in 3rd dimension.
!!  ffti3_local(n3)=local index for 3d dimension
!!  datar(cplex*nfft,nspden)= array in real space.
!!  [action]
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function nctk_write_datar(varname,path,ngfft,cplex,nfft,nspden,&
   comm_fft,fftn3_distrib,ffti3_local,datar,action) result(ncerr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nspden,comm_fft
 character(len=*),intent(in) :: path,varname
 character(len=*),optional,intent(in) :: action
!arrays
 integer,intent(in) :: ngfft(18),fftn3_distrib(ngfft(3)),ffti3_local(ngfft(3))
 real(dp),target,intent(in) :: datar(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: ncid,varid,i3,nproc_fft,me_fft,i3_glob,n1,n2,n3,ispden
 logical :: ionode
 character(len=nctk_slen) :: cplex_name,my_action
 !character(len=500) :: msg
!arrays
 real(dp),allocatable :: glob_datar(:,:)

! *************************************************************************

 ! FIXME: Default should be open but this enters into conflict with the abi_estf stuff!
 ! if we are using MPI-IO since file is not open with HDF5.
 !my_action = "create"; if (present(action)) my_action = action
 my_action = "open"; if (present(action)) my_action = action

 nproc_fft = xmpi_comm_size(comm_fft); me_fft = xmpi_comm_rank(comm_fft)
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)

 ! TODO: Be careful here because we should always create with HDF5 if available
 ! to avoid problems if we have to reread with nproc_fft > 1 and MPI-IO
 ionode = .True.; ncerr = nf90_noerr
 if (nproc_fft == 1) then
   select case(my_action)
   case ("open")
     ncerr = nf90_open(path, mode=nf90_write, ncid=ncid)
   case ("create")
     !ncerr = nf90_create(path, cmode=nf90_clobber, ncid=ncid)
     ncerr = nctk_open_create(ncid, path, comm_fft)
   case default
     MSG_ERROR(sjoin("Wrong action: ", my_action))
   end select

 else
   if (nctk_has_mpiio) then
     call wrtout(std_out, strcat("nctk_write_datar: using MPI-IO to write ", varname, path), "COLL")

     ncerr = nf90_einval
#ifdef HAVE_NETCDF_MPI
     select case(my_action)
     case ("open")
       ncerr = nf90_open(path, mode=nf90_write, comm=comm_fft, info=xmpio_info, ncid=ncid)
     case ("create")
       ncerr = nf90_create(path, cmode=ior(ior(nf90_netcdf4, nf90_mpiio), nf90_write), &
         comm=comm_fft, info=xmpio_info, ncid=ncid)
     case default
       MSG_ERROR(strcat("Wrong action:", my_action))
     end select
#endif
   else
     ! MPI-FFT without MPI-support. Only master does IO
     ionode = (me_fft == master)
     if (ionode) then
       select case(my_action)
       case ("open")
         ncerr = nf90_open(path, mode=nf90_write, ncid=ncid)
       case ("create")
         ncerr = nf90_create(path, cmode=nf90_clobber, ncid=ncid)
       case default
         MSG_ERROR(strcat("Wrong action:", my_action))
       end select
     end if
   end if
 end if
 NCF_CHECK_MSG(ncerr, sjoin("opening file:", path))

 if (ionode) then
   ! Define dims and variables
   !write(std_out,*)"defing dims",trim(varname)," in file: ",path
   cplex_name = strcat("real_or_complex_", varname)
   ncerr = nctk_def_dims(ncid, [&
     nctkdim_t(cplex_name, cplex),&
     nctkdim_t("number_of_grid_points_vector1", n1),&
     nctkdim_t("number_of_grid_points_vector2", n2),&
     nctkdim_t("number_of_grid_points_vector3", n3),&
     nctkdim_t("number_of_components", nspden)], defmode=.True.)
   NCF_CHECK(ncerr)

   ncerr = nctk_def_one_array(ncid, nctkarr_t(name=varname, dtype="dp", shape_str=strcat(cplex_name, &
", number_of_grid_points_vector1, number_of_grid_points_vector2, number_of_grid_points_vector3, number_of_components")),&
   varid=varid)

   ! Add attributes
   varid = nctk_idname(ncid, varname)
   NCF_CHECK(nf90_put_att(ncid, varid, "units", "atomic units"))
   NCF_CHECK(nf90_put_att(ncid, varid, "scale_to_atomic_units", one))
 end if

 if (nproc_fft == 1) then
   ! no MPI-FFT --> write data directly.
   varid = nctk_idname(ncid, varname)
   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, varid, datar, start=[1,1,1,1,1], count=[cplex, n1, n2, n3, nspden]))
   NCF_CHECK(nf90_close(ncid))

 else
   ! Must handle data distribution.
   ABI_CHECK(mod(n3, nproc_fft) == 0, "assuming mod(n3, nproc_fft) == 0")

   i3_glob = -1
   do i3=1,ngfft(3)
     if (fftn3_distrib(i3) == me_fft) then
        i3_glob = i3
        exit
     end if
   end do
   ABI_CHECK(i3_glob > 0, "negative i3_glob")

   !do i3=i3_glob,ngfft(3)
   !  if (fftn3_distrib(i3) /= me_fft) exit
   !  !assert all(ffn3_distrib(i3_glob:i3_glob -1 + ngfft(3) / nproc_fft) == me_fft)
   !end do
   !i3_glob = i3- 1
   !print*,"i3_glob",i3_glob

   if (ionode) then
     NCF_CHECK(nf90_enddef(ncid))
   end if

   ! Array on disk has shape [cplex, n1, n2, n3, nspden]
   if (nctk_has_mpiio) then
     ! Use collective IO.
     ncerr = nf90_einval
     NCF_CHECK(nctk_set_collective(ncid, varid))

     do ispden=1,nspden
       ncerr = nf90_put_var(ncid, varid, datar(:,ispden), start=[1,1,1,i3_glob,ispden], &
                count=[cplex,n1,n2,n3/nproc_fft,1])
       NCF_CHECK(ncerr)
     end do
   else
     ! MPI-FFT without MPI-IO. Collect data (requires more memory and communication)
     ABI_MALLOC(glob_datar, (cplex*product(ngfft(1:3)), nspden))
     call collect_datar(ngfft,cplex,nfft,nspden,datar,comm_fft,fftn3_distrib,ffti3_local,glob_datar,master=master)

     if (ionode) then
       ! Write global array.
       NCF_CHECK(nf90_put_var(ncid, varid, glob_datar, start=[1,1,1,1,1], count=[cplex,n1,n2,n3,nspden]))
     end if
     ABI_FREE(glob_datar)
   end if

   if (ionode) then
     NCF_CHECK(nf90_close(ncid))
   end if
 end if

 !ok = .True.
 ! Sequential IO
 !do rank=0,nproc_fft-1
 !  if (rank == me_fft) then
 !     ncerr = nf90_open(path, mode=nf90_write, ncid=ncid)
 !     do ispden=1,nspden
 !       ncerr = nf90_put_var(ncid, varid, datar(1:,ispden), start=[1,1,1,i3_glob,ispden], &
 !                count=[cplex,ngfft(1),ngfft(2),ngfft(3)/nproc_fft,1])
 !     end do
 !     ncerr = nf90_close(ncid)
 !  end if
 !  call xmpi_barrier(comm_fft)
 !end do

end function nctk_write_datar
!!***

!!****f* m_nctk/nctk_read_datar
!! NAME
!! nctk_read_datar
!!
!! FUNCTION
!!  Read an array in real space in netcdf format
!!
!! INPUTS
!!  path=Filename
!!  varname=Name of the variable to read.
!!  ngfft(18)=information about 3D FFT
!!  cplex=1 for real arrays (e.g. GS rhor), 2 for complex array.
!!  nfft=number of points in the real space FFT mesh treated by this MPI proc
!!  nspden=number of spin-density components
!!  comm_fft=MPI communicator (used only if MPI-FFT).
!!  fftn3_distrib(n3)=rank of the processors which own fft planes in 3rd dimension.
!!  ffti3_local(n3)=local index for 3d dimension
!!  datar(cplex*nfft,nspden)= array in real space.
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function nctk_read_datar(path,varname,ngfft,cplex,nfft,nspden,&
   comm_fft,fftn3_distrib,ffti3_local,datar) result(ncerr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nspden,comm_fft
 character(len=*),intent(in) :: path,varname
!arrays
 integer,intent(in) :: ngfft(18),fftn3_distrib(ngfft(3)),ffti3_local(ngfft(3))
 real(dp),intent(out) :: datar(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: ncid,varid,i3,nproc_fft,me_fft,i3_glob,n1,n2,n3,ispden
 logical :: ionode
!arrays
 real(dp),allocatable :: glob_datar(:,:)

! *************************************************************************

 nproc_fft = xmpi_comm_size(comm_fft); me_fft = xmpi_comm_rank(comm_fft)
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)

 ! TODO: Be careful here because we should always create with HDF5 if available
 ! to avoid problems if we have to reread with nproc_fft > 1 and MPI-IO
 ionode = .True.
 if (nproc_fft == 1) then
   ncerr = nf90_open(path, mode=nf90_nowrite, ncid=ncid)
 else
   if (nctk_has_mpiio) then
     !write(std_out,*)"open_par: ",trim(path)
     !ncerr = nf90_open_par(path, nf90_nowrite,
     ! Don't know why but the format is not autodected!
     ncerr = nf90_einval
#ifdef HAVE_NETCDF_MPI
     ncerr = nf90_open(path, mode=ior(ior(nf90_netcdf4, nf90_mpiio), nf90_nowrite),&
                       comm=comm_fft, info=xmpio_info, ncid=ncid)
#endif
   else
     ! MPI-FFT without MPI-support. Only master does IO
     ionode = (me_fft == master); ncerr = nf90_noerr
     if (ionode) ncerr = nf90_open(path, nf90_nowrite, ncid)
   end if
 end if
 NCF_CHECK_MSG(ncerr, sjoin("opening file: ",path))

 NCF_CHECK(nf90_inq_varid(ncid, varname, varid))
 !write(std_out,*)"about to read varname, ngfft, cplex, nfft, nspden:", trim(varname), ngfft(:3), cplex,nfft,nspden

 ! netcdf array has shape [cplex, n1, n2, n3, nspden]
 if (nproc_fft == 1) then
   ! No MPI-FFT --> easy
   NCF_CHECK(nf90_get_var(ncid, varid, datar, start=[1,1,1,1], count=[cplex, n1, n2, n3, nspden]))
   NCF_CHECK(nf90_close(ncid))

 else
   ! Handle data distribution.
   ABI_CHECK(mod(ngfft(3), nproc_fft) == 0, "assuming mod(n3, nproc_fft) == 0")

   i3_glob = -1
   do i3=1,ngfft(3)
     if (fftn3_distrib(i3) == me_fft) then
        i3_glob = i3
        exit
     end if
   end do
   ABI_CHECK(i3_glob > 0, "negative i3_glob")

   if (nctk_has_mpiio) then
     ! Use parallel IO with collective calls.
     ncerr = nf90_einval
     NCF_CHECK(nctk_set_collective(ncid, varid))

     do ispden=1,nspden
       ncerr = nf90_get_var(ncid, varid, datar(:,ispden), start=[1,1,1,i3_glob,ispden], &
                count=[cplex,ngfft(1),ngfft(2),ngfft(3)/nproc_fft,1])
       NCF_CHECK(ncerr)
     end do
   else
     ! MPI-FFT without MPI-IO. Master read and broadcast (requires more memory and communication)
     ABI_MALLOC(glob_datar, (cplex*product(ngfft(1:3)), nspden))
     if (ionode) then
       NCF_CHECK(nf90_get_var(ncid, varid, glob_datar, start=[1,1,1,1,1], count=[cplex,n1,n2,n3,nspden]))
     end if

     call distrib_datar(ngfft,cplex,nfft,nspden,glob_datar,master,comm_fft,fftn3_distrib,ffti3_local,datar)
     ABI_FREE(glob_datar)
   end if

   if (ionode) then
     NCF_CHECK(nf90_close(ncid))
   end if
 end if

end function nctk_read_datar
!!***

!----------------------------------------------------------------------

!!****f* m_nctk/collect_datar
!! NAME
!!  collect_datar
!!
!! FUNCTION
!! Collect a real-space MPI-FFT distributed array on each proc.
!!
!! INPUTS
!!  ngfft(18)=contain all needed information about 3D FFT (see NOTES at beginning of scfcv)
!!  cplex=1 if real array, 2 for complex
!!  nfft=Number of FFT points treated by this MPI proc
!!  nspden=Second dimension of rhor
!!  rhor(cplex*nfft,nspden)=Array in real space (MPI-FFT distributed)
!!  fftn3_distrib(n3)=rank of the processors which own fft planes in 3rd dimension.
!!  fftn3_local(n3)=local i3 indices
!!  comm_fft=MPI-FFT communicator
!!  [master]=MPI rank, Optional. If present, the global array is available only on master node.
!!
!! OUTPUT
!!   rhor_glob(cplex*nfft_tot,nspden)=Global array
!!
!! PARENTS
!!      m_nctk
!!
!! CHILDREN
!!
!! SOURCE

subroutine collect_datar(ngfft,cplex,nfft,nspden,rhor,comm_fft,fftn3_distrib,ffti3_local,rhor_glob,master)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nspden,comm_fft
 integer,optional,intent(in) :: master
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(in) :: fftn3_distrib(ngfft(3)),ffti3_local(ngfft(3))
 real(dp),intent(in) :: rhor(cplex*nfft,nspden)
 real(dp),intent(out) :: rhor_glob(cplex*product(ngfft(1:3)),nspden)

!Local variables-------------------------------
 integer :: ispden,i1,i2,i3,me_fft,i3_local,my_fftbase,glob_fftbase
 integer :: n1,n2,n3,ierr,nfft_tot

! *************************************************************************

 nfft_tot = product(ngfft(1:3)); me_fft = xmpi_comm_rank(comm_fft)
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)

 if (nfft_tot == nfft) then
   ! full rhor on each node, just do a copy
   rhor_glob = rhor
 else
   ! if MPI-FFT we have to gather the global array on each node.
   rhor_glob = zero
   do ispden=1,nspden
     do i3=1,n3
       if (me_fft == fftn3_distrib(i3)) then
         i3_local = ffti3_local(i3)
         do i2=1,n2
           my_fftbase =   cplex * ( (i2-1)*n1 + (i3_local-1)*n1*n2 )
           glob_fftbase = cplex * ( (i2-1)*n1 + (i3-1)*n1*n2 )
           do i1=1,cplex * n1
             rhor_glob(i1+glob_fftbase,ispden) = rhor(i1+my_fftbase,ispden)
           end do
         end do
       end if
     end do
   end do
   if (present(master)) then
     call xmpi_sum_master(rhor_glob,master,comm_fft,ierr)
   else
     call xmpi_sum(rhor_glob,comm_fft,ierr)
   end if
 end if

end subroutine collect_datar
!!***

!----------------------------------------------------------------------

!!****f* m_nctk/distrib_datar
!! NAME
!!  distrib_datar
!!
!! FUNCTION
!! distribute a real-space MPI-FFT
!!
!! INPUTS
!!  ngfft(18)=contain all needed information about 3D FFT (see NOTES at beginning of scfcv)
!!  cplex=1 if real array, 2 for complex
!!  nfft=Number of FFT points treated by this MPI proc
!!  nspden=Second dimension of rhor
!!  rhor_glob(cplex*nfft_tot,nspden)=Global array
!!  master=The rank of the node that owns the global array.
!!  comm_fft=MPI-FFT communicator
!!  fftn3_distrib(n3)=rank of the processors which own fft planes in 3rd dimension.
!!  fftn3_local(n3)=local i3 indices
!!
!! OUTPUT
!!  rhor(cplex*nfft,nspden)=Array in real space (MPI-FFT distributed)
!!
!! PARENTS
!!      m_nctk
!!
!! CHILDREN
!!
!! SOURCE

subroutine distrib_datar(ngfft,cplex,nfft,nspden,rhor_glob,master,comm_fft,fftn3_distrib,ffti3_local,rhor)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nspden,comm_fft,master
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(in) :: fftn3_distrib(ngfft(3)),ffti3_local(ngfft(3))
 real(dp),intent(out) :: rhor(cplex*nfft,nspden)
 real(dp),intent(inout) :: rhor_glob(cplex*product(ngfft(1:3)),nspden)

!Local variables-------------------------------
 integer :: ispden,i1,i2,i3,me_fft,i3_local,my_fftbase,glob_fftbase
 integer :: n1,n2,n3,ierr,nfft_tot

! *************************************************************************

 nfft_tot = product(ngfft(1:3)); me_fft = xmpi_comm_rank(comm_fft)
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)

 if (nfft_tot == nfft) then
   ! full rhor on each node, just do a copy
   rhor = rhor_glob
 else
   ! if MPI-FFT we have to gather the global array on each node.
   call xmpi_bcast(rhor_glob,master,comm_fft,ierr)
   do ispden=1,nspden
     do i3=1,n3
       if (me_fft == fftn3_distrib(i3)) then
         i3_local = ffti3_local(i3)
         do i2=1,n2
           my_fftbase =   cplex * ( (i2-1)*n1 + (i3_local-1)*n1*n2 )
           glob_fftbase = cplex * ( (i2-1)*n1 + (i3-1)*n1*n2 )
           do i1=1,cplex * n1
             rhor(i1+my_fftbase,ispden) = rhor_glob(i1+glob_fftbase,ispden)
           end do
         end do
       end if
     end do
   end do
 end if

end subroutine distrib_datar
!!***

!----------------------------------------------------------------------

!!****f* m_nctk/var_from_id
!! NAME
!!  var_from_id
!!
!! FUNCTION
!!  Initialize a nctkvar_t object from the variable id
!!
!! INPUTS
!!  ncid=NC file handle
!!  varid=Variable ID
!!
!! OUTPUT
!!  var<nctkvar_t>=Info on the variable.
!!
!! PARENTS
!!      m_nctk
!!
!! CHILDREN
!!
!! SOURCE

subroutine var_from_id(ncid, varid, var)

!Arguments ------------------------------------
 integer, intent(in) :: ncid, varid
 type(nctkvar_t), intent(out) :: var

!Local variables-------------------------------
!scalars
 integer :: ii, ncerr
 !character(len=NF90_MAX_NAME) :: ncname

! *********************************************************************

 ! Get info about the variable.
 var%id = varid
 ncerr = nf90_inquire_variable(ncid, var%id, &
   name=var%name, xtype=var%xtype, ndims=var%ndims, dimids=var%dimids, natts=var%natts)
 NCF_CHECK(ncerr)

 ! Get info about dimensions.
 if (var%ndims > 0) then
   do ii=1,var%ndims
     ncerr = nf90_inquire_dimension(ncid, var%dimids(ii), len=var%dimlens(ii), name=var%dimnames(ii))
     NCF_CHECK(ncerr)
   end do
 end if

 ! Get the number of attributes and their names.
 !if (var%natts > 0) then
 !   do ii=1,var%natts
 !      ncerr = nf90_inq_attname(ncid, var%id, ii, var%attnames(ii))
 !   end do
 !end if

end subroutine var_from_id
!!***

!----------------------------------------------------------------------

!!****f* m_nctk/var_from_name
!! NAME
!!  var_from_name
!!
!! FUNCTION
!!  Initialize a nctkvar_t object from the variable name
!!
!! INPUTS
!!  ncid=NC file handle
!!  name=Variable name
!!
!! OUTPUT
!!  var<nctkvar_t>=Info on the variable.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine var_from_name(ncid, name, var)

!Arguments ------------------------------------
 integer, intent(in) :: ncid
 character(len=*),intent(in) :: name
 type(nctkvar_t), intent(out) :: var

!Local variables-------------------------------
!scalars
 integer :: varid

! *********************************************************************

 varid = nctk_idname(ncid, name)
 call var_from_id(ncid, varid, var)

end subroutine var_from_name
!!***

!!****f* m_nctk/nctk_defwrite_nonana_terms
!! NAME
!! nctk_defwrite_nonana_terms
!!
!! FUNCTION
!!  Write phonon frequencies and displacements for q-->0 in the presence of non-analytical behaviour.
!!
!! INPUTS
!!  ncid=netcdf file id.
!!  iphl2=Index of the q-point to be written to file
!!  nph2l=Number of qpoints.
!!  qph2l(3,nph2l)=List of phonon wavevector directions along which the non-analytical correction
!!    to the Gamma-point phonon frequencies will be calculated
!!    The direction is in CARTESIAN COORDINATES
!!  natom=Number of atoms
!!  phfrq(3*natom)=Phonon frequencies in Ha
!!  cart_displ(2,3*natom,3*natom)=displacements in CARTESIAN coordinates.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      anaddb,m_ifc
!!
!! CHILDREN
!!
!! SOURCE

subroutine nctk_defwrite_nonana_terms(ncid, iphl2, nph2l, qph2l, natom, phfrq, cart_displ, mode)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid,iphl2,nph2l,natom
 character(len=*),intent(in) :: mode
!arrays
 real(dp),intent(in) :: qph2l(3, nph2l)
 real(dp),intent(in) :: phfrq(3*natom)
 real(dp),intent(in) :: cart_displ(2,3*natom,3*natom)

!Local variables-------------------------------
!scalars
 integer :: ncerr, na_phmodes_varid, na_phdispl_varid

! *************************************************************************

 select case (mode)
 case ("define")
   !NCF_CHECK(nctk_def_basedims(ncid, defmode=.True.))
   ncerr = nctk_def_dims(ncid, [nctkdim_t("number_of_non_analytical_directions", nph2l)], defmode=.True.)
   NCF_CHECK(ncerr)

   ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t('non_analytical_directions', "dp", "number_of_cartesian_directions, number_of_non_analytical_directions"),&
   nctkarr_t('non_analytical_phonon_modes', "dp", "number_of_phonon_modes, number_of_non_analytical_directions"),&
   nctkarr_t('non_analytical_phdispl_cart', "dp", &
   "two, number_of_phonon_modes, number_of_phonon_modes, number_of_non_analytical_directions")])
   NCF_CHECK(ncerr)

   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "non_analytical_directions"), qph2l))

 case ("write")

   NCF_CHECK(nf90_inq_varid(ncid, "non_analytical_phonon_modes", na_phmodes_varid))
   NCF_CHECK(nf90_put_var(ncid,na_phmodes_varid,phfrq*Ha_eV,start=[1, iphl2], count=[3*natom, 1]))
   NCF_CHECK(nf90_inq_varid(ncid, "non_analytical_phdispl_cart", na_phdispl_varid))
   ncerr = nf90_put_var(ncid,na_phdispl_varid,cart_displ*Bohr_Ang,&
   start=[1,1,1,iphl2], count=[2,3*natom,3*natom, 1])
   NCF_CHECK(ncerr)

 case default
   MSG_ERROR(sjoin("Wrong value for mode", mode))
 end select

end subroutine nctk_defwrite_nonana_terms
!!***

!!****f* m_nctk/nctk_defwrite_nonana_raman_terms
!! NAME
!! nctk_defwrite_nonana_raman_terms
!!
!! FUNCTION
!! Write the Raman susceptiblities for q-->0 along different directions in the netcdf file.
!!
!! INPUTS
!!  ncid=netcdf file id.
!!  iphl2=Index of the q-point to be written to file.
!!  nph2l=Number of qpoints.
!!  rsus(3*natom,3,3)=List of Raman susceptibilities along the direction corresponding to iphl2.
!!  natom=Number of atoms
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine nctk_defwrite_nonana_raman_terms(ncid, iphl2, nph2l, natom, rsus, mode)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid,natom,iphl2,nph2l
 character(len=*),intent(in) :: mode
!arrays
 real(dp),intent(in) :: rsus(3*natom,3,3)

!Local variables-------------------------------
!scalars
 integer :: ncerr, raman_sus_varid

! *************************************************************************

!Fake use of nph2l, to keep it as argument. This should be removed when nph2l will be used.
 if(.false.)then
  ncerr=nph2l
 endif

 select case (mode)
 case ("define")
   NCF_CHECK(nctk_def_basedims(ncid, defmode=.True.))
   ncerr = nctk_def_arrays(ncid, [ nctkarr_t("non_analytical_raman_sus", "dp", &
"number_of_non_analytical_directions,number_of_phonon_modes,number_of_cartesian_directions,number_of_cartesian_directions")])
   NCF_CHECK(ncerr)

   NCF_CHECK(nctk_set_datamode(ncid))

 case ("write")

   NCF_CHECK(nf90_inq_varid(ncid, "non_analytical_raman_sus", raman_sus_varid))
   ncerr = nf90_put_var(ncid,raman_sus_varid,rsus,&
     start=[iphl2,1,1,1], count=[1,3*natom,3,3])
   NCF_CHECK(ncerr)

 case default
   MSG_ERROR(sjoin("Wrong value for mode", mode))
 end select

end subroutine nctk_defwrite_nonana_raman_terms
!!***

!!****f* m_nctk/nctk_defwrite_raman_terms
!! NAME
!! nctk_defwrite_raman_terms
!!
!! FUNCTION
!! Write the Raman susceptiblities for q=0 and also the phonon frequncies at gamma.
!!
!! INPUTS
!!  ncid=netcdf file id.
!!  rsus(3*natom,3,3)=List of Raman susceptibilities.
!!  natom=Number of atoms
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine nctk_defwrite_raman_terms(ncid, natom, rsus, phfrq)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid,natom
!arrays
 real(dp),intent(in) :: rsus(3*natom,3,3)
 real(dp),intent(in) :: phfrq(3*natom)

!Local variables-------------------------------
!scalars
 integer :: ncerr, raman_sus_varid, phmodes_varid

! *************************************************************************

 NCF_CHECK(nctk_def_basedims(ncid, defmode=.True.))
 ncerr = nctk_def_arrays(ncid, [ nctkarr_t("raman_sus", "dp", &
  "number_of_phonon_modes,number_of_cartesian_directions,number_of_cartesian_directions"), &
  nctkarr_t("gamma_phonon_modes", "dp", "number_of_phonon_modes")])
 NCF_CHECK(ncerr)

 NCF_CHECK(nctk_set_datamode(ncid))

 NCF_CHECK(nf90_inq_varid(ncid, "raman_sus", raman_sus_varid))
 NCF_CHECK(nf90_put_var(ncid,raman_sus_varid,rsus))
 NCF_CHECK(nf90_inq_varid(ncid, "gamma_phonon_modes", phmodes_varid))
 NCF_CHECK(nf90_put_var(ncid,phmodes_varid,phfrq*Ha_eV))

end subroutine nctk_defwrite_raman_terms
!!***

#endif

!!****f* m_nctk/create_nc_file
!! NAME
!! create_nc_file
!!
!! FUNCTION
!! Create an NetCDF file including a dimension one definition
!!
!! INPUTS
!!
!! OUTPUT
!!
!! TODO:
!!  Remove
!!
!! PARENTS
!!      outvars
!!
!! CHILDREN
!!
!! SOURCE

subroutine create_nc_file (filename,ncid)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: ncid
!arrays
character(len=*),intent(in) :: filename

!Local variables-------------------------------
#if defined HAVE_NETCDF
integer :: one_id
integer :: ncerr
#endif

! *************************************************************************

 ncid = 0
#if defined HAVE_NETCDF
!Create the NetCDF file
 ncerr=nf90_create(path=filename,cmode=NF90_CLOBBER,ncid=ncid)
 NCF_CHECK_MSG(ncerr, sjoin('Error while creating:', filename))
 ncerr=nf90_def_dim(ncid,'one',1,one_id)
 NCF_CHECK_MSG(ncerr,'nf90_def_dim')
#endif

 end subroutine create_nc_file
!!***

!!****f* m_nctk/write_var_netcdf
!!
!! NAME
!! write_var_netcdf
!!
!! FUNCTION
!! Write variable into a netcdf dataset
!!
!! TODO:
!!  Remove!
!!
!! INPUTS
!! arr_int
!! arr_real
!! marr
!! narr
!! typevar
!! varname
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      prttagm,prttagm_images
!!
!! CHILDREN
!!
!! SOURCE

subroutine write_var_netcdf(arr_int,arr_real,marr,narr,ncid,typevar,varname)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: narr,marr,ncid
 character(len=*),intent(in) :: varname
 character(len=3),intent(in) :: typevar
!arrays
 integer,intent(in) :: arr_int(marr)
 real(dp),intent(in) :: arr_real(marr)

!Local variables-------------------------------
!scalars
 integer :: var_id,var_type,vardim_id,ncerr
 !character(len=500) :: msg

! *************************************************************************

 !write(std_out,*)"about to write varname: ",trim(varname)

#if defined HAVE_NETCDF
 if (ncid>0) then
!  ### Put the file in definition mode
   ncerr=nf90_redef(ncid)
   if (ncerr/=NF90_NOERR.and.ncerr/=NF90_EINDEFINE) then
     NCF_CHECK_MSG(ncerr,'nf90_redef')
   end if
!  ### Define the dimensions
   if (narr==1)then
     ncerr=nf90_inq_dimid(ncid,'one',vardim_id)
     NCF_CHECK_MSG(ncerr,'nf90_inq_varid')
   else
     ncerr=nf90_def_dim(ncid,trim(varname),narr,vardim_id)
     NCF_CHECK_MSG(ncerr,'nf90_def_dim')
   end if
!  ### Define the variables
   if (typevar=='INT') then
     var_type=NF90_INT
   else if (typevar=='DPR') then
     var_type=NF90_DOUBLE
   end if
   ncerr=nf90_def_var(ncid, trim(varname), var_type, vardim_id, var_id)
   NCF_CHECK_MSG(ncerr,'nf90_def_var')
!  ### Put the file in data mode
   ncerr=nf90_enddef(ncid)
   if (ncerr/=NF90_NOERR.and.ncerr/=NF90_ENOTINDEFINE) then
     NCF_CHECK_MSG(ncerr,'nf90_enddef')
   end if
!  ### Write variables into the dataset
   if (typevar=='INT') then
     ncerr=nf90_put_var(ncid,var_id,arr_int,start=(/1/),count=(/narr/))
   else if (typevar=='DPR') then
     ncerr=nf90_put_var(ncid,var_id,arr_real,start=(/1/),count=(/narr/))
   end if
   NCF_CHECK_MSG(ncerr,'nf90_put_var')
 end if
#endif

end subroutine write_var_netcdf
!!***

!!****f* ABINIT/write_eig
!!
!! NAME
!! write_eig
!!
!! FUNCTION
!! Write the eigenvalues band by band and k point by k point
!! in a NetCDF file format
!!
!! INPUTS
!! filname = Filename of the file where the history will be stored
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      clnup1
!!
!! CHILDREN
!!      ab_define_var
!!
!! SOURCE

subroutine write_eig(eigen,filename,kptns,mband,nband,nkpt,nsppol)

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: filename
 integer,intent(in) :: nkpt,nsppol,mband
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
 real(dp),intent(in) :: kptns(3,nkpt)

!Local variables-------------------------------
!scalars
 integer :: ncerr,ncid,ii
 integer :: xyz_id,nkpt_id,mband_id,nsppol_id
 integer :: eig_id,kpt_id,nbk_id,nbk
 integer :: ikpt,isppol,nband_k,band_index
 real(dp):: convrt
!arrays
 integer :: dimEIG(3),dimKPT(2),dimNBK(2)
 integer :: count2(2),start2(2)
 integer :: count3(3),start3(3)
 real(dp):: band(mband)

! *********************************************************************

#if defined HAVE_NETCDF

 convrt=1.0_dp

!1. Create netCDF file
 ncerr = nf90_create(path=trim(filename),cmode=NF90_CLOBBER, ncid=ncid)
 NCF_CHECK_MSG(ncerr," create netcdf EIG file")

!2. Define dimensions
 ncerr = nf90_def_dim(ncid,"xyz",3,xyz_id)
 NCF_CHECK_MSG(ncerr," define dimension xyz")

 ncerr = nf90_def_dim(ncid,"mband",mband,mband_id)
 NCF_CHECK_MSG(ncerr," define dimension mband")

 ncerr = nf90_def_dim(ncid,"nkpt",nkpt,nkpt_id)
 NCF_CHECK_MSG(ncerr," define dimension nkpt")

 ncerr = nf90_def_dim(ncid,"nsppol",nsppol,nsppol_id)
 NCF_CHECK_MSG(ncerr," define dimension nsppol")

!Dimensions for EIGENVALUES
 dimEIG = (/ mband_id, nkpt_id, nsppol_id /)
!Dimensions for kpoint positions
 dimKPT = (/ xyz_id, nkpt_id /)
!Dimensions for number kpoints per band and spin
!dimNBK = (/ nkpt_id, nsppol_id /)
 dimNBK = (/ nkpt_id, nsppol_id /)

!3. Define variables

 call ab_define_var(ncid, dimEIG, eig_id, NF90_DOUBLE,&
& "Eigenvalues",&
& "Values of eigenvalues",&
& "Hartree")
 call ab_define_var(ncid, dimKPT, kpt_id, NF90_DOUBLE,"Kptns",&
& "Positions of K-points in reciprocal space",&
& "Dimensionless")
 call ab_define_var(ncid, dimNBK, nbk_id, NF90_INT,"NBandK",&
& "Number of bands per kpoint and Spin",&
& "Dimensionless")

!4. End define mode
 ncerr = nf90_enddef(ncid)
 NCF_CHECK_MSG(ncerr," end define mode")

!5 Write kpoint positions
 do ikpt=1,nkpt
   start2 = (/ 1, ikpt /)
   count2 = (/ 3, 1 /)
   ncerr = nf90_put_var(ncid, kpt_id,&
&   kptns(1:3,ikpt),&
&   start = start2,&
&   count = count2)
   NCF_CHECK_MSG(ncerr," write variable kptns")
 end do


!6 Write eigenvalues
 band_index=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     start3 = (/ 1, ikpt, isppol /)
     count3 = (/ mband, 1, 1 /)
     band(:)=zero
     do ii=1,nband_k
       band(ii)=eigen(band_index+ii)
     end do
     ncerr = nf90_put_var(ncid, eig_id,&
&     band,&
&     start = start3,&
&     count = count3)
     NCF_CHECK_MSG(ncerr," write variable band")

     band_index=band_index+nband_k
   end do
 end do

!6 Write Number of bands per kpoint and Spin

 do isppol=1,nsppol
   do ikpt=1,nkpt
     start2 = (/ ikpt, 1 /)
     count2 = (/ 1, 1 /)
     nbk=nband(ikpt+(isppol-1)*nkpt)
     ncerr = nf90_put_var(ncid, nbk_id,&
&     nbk,&
&     start = start2)
     NCF_CHECK_MSG(ncerr," write variable nband")
   end do
 end do

!7 Close file

 ncerr = nf90_close(ncid)
 NCF_CHECK_MSG(ncerr," close netcdf EIG file")
#endif

end subroutine write_eig
!!***

END MODULE m_nctk
!!***
