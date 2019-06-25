!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_io_screening
!! NAME
!!  m_io_screening
!!
!! FUNCTION
!!  This module contains the definition of the header of the
!!  _SCR and _SUSC file as well as methods used to read/write/echo.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2019 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_io_screening

 use defs_basis
 use defs_abitypes
 use m_abicore
#if defined HAVE_MPI2
 use mpi
#endif
 use m_xmpi
 use m_mpiotk
 use m_nctk
 use m_errors
 use iso_c_binding
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_hdr
 use m_sort

 use m_gwdefs,          only : em1params_t, GW_TOLQ
 use m_fstrings,        only : sjoin, itoa, endswith, replace_ch0
 use m_copy,            only : alloc_copy
 use m_io_tools,        only : open_file, file_exists, iomode2str
 use m_numeric_tools,   only : print_arr, remove_copies, imax_loc
 use m_bz_mesh,         only : isequalk

 implicit none

 private
!!***

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

 character(len=nctk_slen),public,parameter :: e_ncname="dielectric_function"
 character(len=nctk_slen),public,parameter :: em1_ncname="inverse_dielectric_function"
 character(len=nctk_slen),public,parameter :: chi0_ncname="polarizability"

 public :: ncname_from_id       ! return the name of the netcdf variable from the id

!!****t* m_io_screening/hscr_t
!! NAME
!!  hscr_t
!!
!! FUNCTION
!!  The structure defining the header of the SCR/SUSC file.
!!  hscr_t contains the most important dimensions associated to the SCR/SUSC matrix,
!!  important GW metadata and the Abint header. The SCR/SUS matrices are saved
!!  with the same format. There are nqibz blocks, each block contains (npwe,npwe,nomega) matrices.
!!  SCR and SUS files mainly differ for what concerns the treatment of the q-->0 limit.
!!  The treatment of the non-analytic behaviour is not yet implemented but the main ideas are
!!  sketched in the NOTES below.
!!
!! NOTES
!!  On the treatment of the q-->0 limit
!!
!!  1) q=Gamma should be the first q-point
!!
!!  2) This point contains an initial section with data used to treat the q-->0 limit, followed
!!     by the SCR/SUS matrix evaluated for the small q-point (qlwl). The header should contains
!!     enough info so that we can skip this section and use the pre-existing routines to read
!!     the matrices.
!!
!!  3) The data stored in the q-->0 section depends on the type of matrix stored in the file:
!!
!!     SUS file: we store the tensor and the wings needed to reconstruct the q-dependence around Gamma.
!!       This section if followed by the X(G,G') matrix evaluated at qlwl.
!!       Note that the content of the SUS file does not depend on a possible cutoff in vcoul.
!!
!!     SCR file: two cases must be considered:
!!
!!     No cutoff in vcoul:
!!        The q-->0 section stores the tensor, the wings as well as the inverse of the body B^{-1}.
!!        This info is used to integrate W(q) around q==Gamma.
!!
!!     cutoff in vcoul:
!!        The content of the q-->0 section depends on dimensionality of the system.
!!        In 2D we need info on chi0 as well as A and a. See http://arxiv.org/pdf/1511.00129v1.pdf
!!        I think that this kind of calculations are easy to implement if we start from the SUS file.
!!        Ok, we have to recompute e-1 at each run but the logic is easier to implement.
!!
!! SOURCE

 type,public :: hscr_t

  integer :: id
    ! Matrix identifier: 1 for chi0, 2 for chi, 3 for epsilon, 4 for espilon^{-1}

  integer :: ikxc
    ! Kxc kernel used,
    ! 0 for None (RPA), >0 for static TDDFT (=ixc), <0 for frequency-dependent TDDFT

  integer :: inclvkb
    ! q-->0 treatment, 0 for None, 1-2 for transversal gauge, 3 for longitudinal

  integer :: headform
    ! format of the SCR header

  integer :: fform
    ! File format

  integer :: gwcalctyp
    ! Calculation type (G0W0, G0W, GW ...)

  integer :: nI,nJ
    ! Number of spin components (rows,columns) in chi|eps^-1. (1,1) if collinear.
    ! The internal representation of the matrix is eps(nI*npwe,nJ*npwe)

  integer :: nqibz
    ! Number of q-points in the IBZ.

  integer :: nqlwl
    ! Number of points for the treatment of the long wavelength limit.

  integer :: nomega
    ! Total number of frequencies.

  integer :: nbnds_used
    ! Number of bands used during the screening calculation (only for info)

  integer :: npwe
    ! Number of G vectors reported on the file.

  integer :: npwwfn_used
    ! Number of G vectors for wavefunctions used during the screening calculation (only for info)

  integer :: spmeth
    ! Method used to approximate the delta function in the expression for Im Chi_0

  integer :: test_type
    ! 1 for TEST-PARTICLE, 2 for TEST-ELECTRON.

  integer :: tordering
    ! 1 for Time-Ordered, 2 for Advanced, 3 for Retarded.

! HSCR_NEW
  integer :: awtr
  ! Input variable (time-reversal symmetry in RPA expression)

  integer :: icutcoul
  ! Input variable (Coulomb cutoff)

  integer :: gwcomp
  ! Input variable (GW compensation energy technique)

  integer :: gwgamma
  ! Input variable Vertex correction
! HSCR_NEW

  real(dp) :: mbpt_sciss
    ! Scissor Energy, zero if not used

  real(dp) :: spsmear
    ! Smearing of the delta in case of spmeth==2

  real(dp) :: zcut
    ! Imaginary shift to avoid the poles along the real axis.

! HSCR_NEW
  real(dp) :: gwencomp
   ! Input variable (GW compensation energy technique)

  character(len=3) :: kind_cdata
  ! Flag to signal whether the data is in single or double precision ("spc" or "dpc")
  ! For the time being, we always write/read in double precision.
  ! This flag could be use to reduce the memory requirements if spc:
  ! we run calculations in single precision dump the results with the same precision without
  ! having to allocate extra memory.
! HSCR_NEW

!arrays

! HSCR_NEW
  real(dp) :: vcutgeo(3)
   ! Input variable (defines coulomb cutoff)
! HSCR_NEW

  character(len=80) :: titles(2)
    ! Titles describing the content of the file.

  integer,allocatable  :: gvec(:,:)
    ! gvec(3,npwe)
    ! G vectors in reduced coordinates.

  real(dp),allocatable :: qibz(:,:)
    ! qibz(3,nqibz)
    ! q-points in the IBZ in reduced coordinates.

  real(dp),allocatable :: qlwl(:,:)
    ! qlwl(3,nqlwl)
    ! q-points for the long wave-length limit treatment (r.l.u)

  complex(dpc),allocatable :: omega(:)
    ! omega(nomega)
    ! All frequencies calculated both along the real and the imaginary axis.
    ! Real frequencies are packed in the first section.
    ! TODO: Add frequency mesh type?

  type(hdr_type) :: hdr
    ! The abinit header.

 end type hscr_t
!!***

 integer,private,parameter :: HSCR_KNOWN_HEADFORMS(1) = [80]
 ! The list of headforms used for SCR/SUSC so far.

 integer,public,parameter :: HSCR_LATEST_HEADFORM = HSCR_KNOWN_HEADFORMS(size(HSCR_KNOWN_HEADFORMS))
 ! The latest headform used when writing.

 public :: hscr_from_file       ! Read the header from file.
 public :: hscr_io              ! I/O of the header (read/write/echo).
 !public :: hscr_fort_read
 !public :: hscr_fort_write
 !public :: hscr_ncwread
 !public :: hscr_ncwrite
 !public :: hscr_echo            ! I/O of the header (read/write/echo).
 public :: hscr_print           ! Print the SCR related part of the header.
 public :: hscr_new             ! Create header.
 public :: hscr_bcast           ! Transmit the header.
 public :: hscr_free            ! Free the header.
 public :: hscr_copy            ! Copy the SCR|SUSC header.
 public :: hscr_merge           ! Merge two or more headers.
 public :: write_screening      ! Write a q-slice of the matrix in G-space.
 public :: read_screening       ! Read the content of the (SCR|SUSC) file placed after the header.

! Tools used in mrgscr.
 public :: ioscr_qmerge         ! Produce new file by merging the q-points stored in other files.
 public :: ioscr_qrecover       ! Recover q-points from a corrupted file produced e.g. from an interrupted run
 public :: ioscr_wmerge         ! Produce new file by merging the frequencies stored in other files.
 public :: ioscr_wremove        ! Produce new file by removing selected frequencies in the initial file.
 !public :: ioscr_nc2fort       ! Convert a netcdf file to a Fortran file. TODO

CONTAINS  !================================================================================================
!!***

!!****f* m_io_screening/ncname_from_id
!! NAME
!!  ncname_from_id
!!
!! FUNCTION
!!  Return the name of the netcdf variable (chi0, espilon, em1...) from the id.
!!
!! SOURCE

character(len=nctk_slen) function ncname_from_id(id) result(varname)

  integer,intent(in) :: id

  varname = "None"
  if (id == 1) varname = chi0_ncname
  if (id == 3) varname = e_ncname
  if (id == 4) varname = em1_ncname
  ABI_CHECK(varname /= "None", "Wrong id")

end function ncname_from_id
!!***

!!****f* m_io_screening/hscr_from_file
!! NAME
!!  hscr_from_file
!!
!! FUNCTION
!!  Read the header of the (SCR/SUS) file
!!
!! INPUTS
!!  path=File name
!!  comm = MPI communicator.
!!
!! OUTPUT
!!  hscr<hscr_t>=The header.
!!  fform=Kind of the array in the file (0 signals an error)
!!
!! PARENTS
!!      m_io_screening,m_screen,m_screening,mrgscr,setup_bse
!!
!! CHILDREN
!!      hscr_copy,hscr_io,hscr_print,read_screening,write_screening,wrtout
!!
!! SOURCE

subroutine hscr_from_file(hscr,path,fform,comm)

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: path
 integer,intent(in) :: comm
 integer,intent(out) :: fform
 type(hscr_t),intent(out) :: hscr

!Local variables-------------------------------
!scalars
 integer,parameter :: rdwr5=5,master=0
 integer :: unt,my_rank,ierr
 character(len=500) :: msg

! *************************************************************************

 my_rank = xmpi_comm_rank(comm)

 ! Master reads and broadcasts.
 if (my_rank == master) then
   if (.not. endswith(path, ".nc")) then
     ! Fortran-IO
     if (open_file(path,msg,newunit=unt,form="unformatted", status="old",action="read") /= 0) then
       MSG_ERROR(msg)
     end if
     call hscr_io(hscr,fform,rdwr5,unt,xmpi_comm_self,master,IO_MODE_FORTRAN)
     close(unt)
   else
     ! Netcdf format
#ifdef HAVE_NETCDF
     NCF_CHECK(nctk_open_read(unt, path, xmpi_comm_self))
     call hscr_io(hscr,fform,rdwr5,unt,xmpi_comm_self,master,IO_MODE_ETSF)
     NCF_CHECK(nf90_close(unt))
#else
     NETCDF_NOTENABLED_ERROR()
#endif
   end if

   ABI_CHECK(fform /= 0, sjoin("hscr_io returned fform == 0 while reading:", path))
 end if

 ! Broadcast data.
 if (xmpi_comm_size(comm) > 1) then
   call hscr_bcast(hscr,master,my_rank,comm)
   call xmpi_bcast(fform,master,comm,ierr)
 end if

end subroutine hscr_from_file
!!***

!!****f* m_io_screening/hscr_io
!! NAME
!!  hscr_io
!!
!! FUNCTION
!! This subroutine deals with the I/O of the hscr_t structured variables (read/write/echo).
!! According to the value of rdwr, it reads the header of a file, writes it, or echo the value
!! of the structured variable to a file. Note that, when reading, different records of hscr_t
!! are allocated here, according to the values of the read variables. Records of hscr_t should be
!! deallocated correctly by a call to hdr_free when hscr_t is not used anymore.
!!
!! INPUTS
!!  iomode=Option defining the file format of the external file.
!!  comm=MPI communicator.
!!  master=rank of the master node in comm, usually 0
!!  rdwr= if 1, read the hscr_t structured variable from the header of the file,
!!        if 2, write the header to unformatted file
!!        if 3, echo part of the header to formatted file (records 1 and 2)
!!        if 4, echo the header to formatted file
!!        if 5, read the hscr_t without rewinding (unformatted)
!!        if 6, write the hscr_t without rewinding (unformatted)
!!  unt=unit number of the file (unformatted if rdwr=1, 2, 5 or 6 formatted if rdwr=3,4)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  The following variables are both input or output :
!!  fform=kind of the array in the file
!!   if rdwr=1,5 : will be output ; if the reading fail, return fform=0
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!  hscr_t <type(hscr_t)>=the header structured variable
!!   if rdwr=1,5 : will be output
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!
!! NOTES
!! In all cases, the file is supposed to be open already
!! When reading (rdwr=1) or writing (rdwr=2), rewind the file
!! When echoing (rdwr=3) does not rewind the file.
!! When reading (rdwr=5) or writing (rdwr=6), DOES NOT rewind the file
!!
!! In writing mode, the routine is supposed to called by the master node.
!! no check is done, it is up to the developer.
!!
!! PARENTS
!!      m_io_screening,m_screening,screening
!!
!! CHILDREN
!!      hscr_copy,hscr_io,hscr_print,read_screening,write_screening,wrtout
!!
!! SOURCE

subroutine hscr_io(hscr,fform,rdwr,unt,comm,master,iomode)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(inout) :: fform
 integer,intent(in) :: rdwr,unt,iomode,comm,master
 type(hscr_t),target,intent(inout) :: hscr

!Local variables-------------------------------
!scalars
 integer :: my_rank,nprocs,ncerr,ncid,varid,ierr !ii
 character(len=500) :: errmsg
 character(len=nctk_slen) :: varname,head_shape,wing_shape
!arrays
 real(dp),allocatable :: real_omega(:,:)
 real(dp), ABI_CONTIGUOUS pointer :: r2vals(:,:) !,rvals3(:,:,:)
! *************************************************************************

 DBG_ENTER("COLL")
 !@hscr_t
 ABI_UNUSED(master) ! FIXME

 ! Initialize MPI info for comm ===
 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 if (rdwr==1 .or. rdwr==5) then

   if (.True.) then
   ! TODO: only master should read but then I have to skip the header.
   !if (my_rank == master) then
     ! Read the abinit header, rewinding of the file (if any) is done here.
     if (iomode==IO_MODE_FORTRAN) then
       call hdr_fort_read(hscr%hdr, unt, fform, rewind=(rdwr==1))
     else if (iomode==IO_MODE_ETSF) then
       call hdr_ncread(hscr%hdr, unt, fform)
     end if

     ! Reset the variables absent in old versions.
     Hscr%fform=fform

     if (iomode==IO_MODE_FORTRAN .or. iomode==IO_MODE_MPI) then
       select case (fform)
       case (1003, 1004)
         ! File format for epsilon^-1, espilon, chi0
         read(unt, err=10, iomsg=errmsg)hscr%titles
         read(unt, err=10, iomsg=errmsg)&
           hscr%id, hscr%ikxc, hscr%inclvkb, hscr%headform, hscr%fform, hscr%gwcalctyp,&
           hscr%nI, hscr%nJ, hscr%nqibz, hscr%nqlwl, hscr%nomega, hscr%nbnds_used,&
           hscr%npwe, hscr%npwwfn_used, hscr%spmeth, hscr%test_type, hscr%tordering,&
           hscr%awtr, hscr%icutcoul, hscr%gwgamma, hscr%vcutgeo(1:3)

         ! Read real scalars
         read(unt, err=10, iomsg=errmsg)&
           hscr%mbpt_sciss, hscr%spsmear, hscr%zcut, hscr%gwencomp,hscr%kind_cdata

         ! Allocate arrays and read them
         call hscr_malloc(hscr, hscr%npwe, hscr%nqibz, hscr%nomega, hscr%nqlwl)

         read(unt, err=10, iomsg=errmsg)hscr%gvec(:,:)
         read(unt, err=10, iomsg=errmsg)hscr%qibz(:,:)
         read(unt, err=10, iomsg=errmsg)hscr%omega(:)

         ! Read data for q-->0 limit.
         if (hscr%nqlwl>0) then
           read(unt, err=10, iomsg=errmsg)hscr%qlwl(:,:)
         end if

       case default
         MSG_BUG(sjoin('Wrong fform read:', itoa(fform)))
       end select

     else if (iomode==IO_MODE_ETSF) then
#ifdef HAVE_NETCDF
       ncid = unt

       select case (fform)
       case (1003, 1004)
         ! Get dimensions and allocate arrays.
         NCF_CHECK(nctk_get_dim(ncid, "number_of_coefficients_dielectric_function", hscr%npwe))
         NCF_CHECK(nctk_get_dim(ncid, "number_of_qpoints_dielectric_function", hscr%nqibz))
         NCF_CHECK(nctk_get_dim(ncid, "number_of_frequencies_dielectric_function", hscr%nomega))
         NCF_CHECK(nctk_get_dim(ncid, "number_of_qpoints_gamma_limit", hscr%nqlwl))
         NCF_CHECK(nctk_get_dim(ncid, "nI", hscr%ni))
         NCF_CHECK(nctk_get_dim(ncid, "nJ", hscr%nj))
         call hscr_malloc(hscr, hscr%npwe, hscr%nqibz, hscr%nomega, hscr%nqlwl)

         varid = nctk_idname(ncid, 'reduced_coordinates_plane_waves_dielectric_function')
         NCF_CHECK(nf90_get_var(ncid, varid, hscr%gvec, start=[1,1,1]))
         NCF_CHECK(nf90_get_var(ncid, vid('qpoints_dielectric_function'), hscr%qibz))

         ABI_MALLOC(real_omega, (2, hscr%nomega))
         NCF_CHECK(nf90_get_var(ncid, vid('frequencies_dielectric_function'), real_omega))
         hscr%omega = dcmplx(real_omega(1,:), real_omega(2,:))
         ABI_FREE(real_omega)

         ! Read extra data added in new header.
         NCF_CHECK(nf90_get_var(ncid, vid("vcutgeo"), hscr%vcutgeo))
         NCF_CHECK(nf90_get_var(ncid, vid("id"), hscr%id))
         NCF_CHECK(nf90_get_var(ncid, vid("ikxc"), hscr%ikxc))
         NCF_CHECK(nf90_get_var(ncid, vid("inclvkb"), hscr%inclvkb))
         NCF_CHECK(nf90_get_var(ncid, vid("headform"), hscr%headform))
         NCF_CHECK(nf90_get_var(ncid, vid("fform"), hscr%fform))
         NCF_CHECK(nf90_get_var(ncid, vid("gwcalctyp"), hscr%gwcalctyp))
         NCF_CHECK(nf90_get_var(ncid, vid("nbands_used"), hscr%nbnds_used))
         NCF_CHECK(nf90_get_var(ncid, vid("npwwfn_used"), hscr%npwwfn_used))
         NCF_CHECK(nf90_get_var(ncid, vid("spmeth"), hscr%spmeth))
         NCF_CHECK(nf90_get_var(ncid, vid("test_type"), hscr%test_type))
         NCF_CHECK(nf90_get_var(ncid, vid("tordering"), hscr%tordering))
         NCF_CHECK(nf90_get_var(ncid, vid("awtr"), hscr%awtr))
         NCF_CHECK(nf90_get_var(ncid, vid("icutcoul"), hscr%icutcoul))
         NCF_CHECK(nf90_get_var(ncid, vid("gwcomp"), hscr%gwcomp))
         NCF_CHECK(nf90_get_var(ncid, vid("gwgamma"), hscr%gwgamma))
         NCF_CHECK(nf90_get_var(ncid, vid("mbpt_sciss"), hscr%mbpt_sciss))
         NCF_CHECK(nf90_get_var(ncid, vid("spsmear"), hscr%spsmear))
         NCF_CHECK(nf90_get_var(ncid, vid("zcut"), hscr%zcut))
         NCF_CHECK(nf90_get_var(ncid, vid("gwencomp"), hscr%gwencomp))
         NCF_CHECK(nf90_get_var(ncid, vid("kind_cdata"), hscr%kind_cdata))
         call replace_ch0(hscr%kind_cdata)

         NCF_CHECK(nf90_get_var(ncid, vid("titles"), hscr%titles))
         call replace_ch0(hscr%titles(:))

         ! TODO Read it
         if (hscr%nqlwl /= 0) then
           NCF_CHECK(nf90_get_var(ncid, vid("qpoints_gamma_limit"), hscr%qlwl))
         end if

       case default
         MSG_BUG(sjoin('Unsupported fform read:',itoa(fform)))
       end select
#endif
     else
       MSG_ERROR(sjoin("Unsupported value of iomode:", iomode2str(iomode)))
     end if

   end if ! master

   !call hscr_bcast(hscr, master, my_rank, comm)
   !call hscr_mpio_skip(mpio_fh,fform,offset)

 else if (rdwr==2.or.rdwr==6) then
   ! Writing the header of an unformatted file.
   ! Always use the latest version.

   if (iomode==IO_MODE_FORTRAN .or. iomode==IO_MODE_MPI) then
     ! Write the abinit header.
     call hdr_fort_write(hscr%hdr, unt, fform, ierr)
     ABI_CHECK(ierr == 0, "hdr_fort_write retured ierr != 0")

     write(unt, err=10, iomsg=errmsg)hscr%titles

     ! Write integers
     write(unt, err=10, iomsg=errmsg)&
       hscr%id, hscr%ikxc, hscr%inclvkb, hscr%headform, hscr%fform, hscr%gwcalctyp,&
       hscr%nI, hscr%nJ, hscr%nqibz, hscr%nqlwl, hscr%nomega, hscr%nbnds_used,&
       hscr%npwe, hscr%npwwfn_used, hscr%spmeth, hscr%test_type, hscr%tordering,&
       hscr%awtr, hscr%icutcoul, hscr%gwgamma, hscr%vcutgeo(1:3)

     ! Write real scalars
     write(unt, err=10, iomsg=errmsg)&
       hscr%mbpt_sciss, hscr%spsmear, hscr%zcut, hscr%gwencomp,hscr%kind_cdata

     ! Write arrays
     write(unt, err=10, iomsg=errmsg)hscr%gvec(:,:)
     write(unt, err=10, iomsg=errmsg)hscr%qibz(:,:)
     write(unt, err=10, iomsg=errmsg)hscr%omega(:)

     ! Add q-points for heads and wings for q-->0.
     if (hscr%nqlwl>0) then
       write(unt, err=10, iomsg=errmsg)hscr%qlwl(:,:)
     end if

   else if (iomode == IO_MODE_ETSF) then
#ifdef HAVE_NETCDF
     ncid = unt
     ! Write the abinit header, rewinding of the file (if any) is done here.
     NCF_CHECK(hdr_ncwrite(hscr%hdr, ncid, fform, nc_define=.True.))

     ! Define dimensions
     ! Part 2) of etsf-io specifications
     ! FIXME: Spin is only used in particular cases, We usuall get the trace of W in spin space
     ! and I'm not gonna allocate extra memory just to have up up, down down
     ! Besides number_of_spins should be replaced by `number_of_spins_dielectric_function`
     ! Should add spin_dependent attribute.
     ncerr = nctk_def_dims(ncid, [&
       nctkdim_t("complex", 2), nctkdim_t("number_of_reduced_dimensions", 3),&
       nctkdim_t("number_of_frequencies_dielectric_function", hscr%nomega), &
       nctkdim_t("number_of_qpoints_dielectric_function", hscr%nqibz),&
       nctkdim_t("number_of_qpoints_gamma_limit", hscr%nqlwl),&
       nctkdim_t("number_of_spins", hscr%hdr%nsppol),&
       nctkdim_t("nI", hscr%nI), nctkdim_t("nJ", hscr%nJ),&
       nctkdim_t("number_of_coefficients_dielectric_function", hscr%npwe)], defmode=.True.)
     NCF_CHECK(ncerr)

     ! Part 3) of the specs
     ! (note that, in the specs, the Gs depend on the q-point but npwe is a scalar
     ! basis_set is added by the abinit header.
     ! FIXME: g-vectors are not written properly.
     ncerr = nctk_def_arrays(ncid, [&
       ! Standard
       nctkarr_t('frequencies_dielectric_function', "dp", 'complex, number_of_frequencies_dielectric_function'), &
       nctkarr_t('qpoints_dielectric_function', "dp", 'number_of_reduced_dimensions, number_of_qpoints_dielectric_function'),&
       nctkarr_t('qpoints_gamma_limit', "dp", 'number_of_reduced_dimensions, number_of_qpoints_gamma_limit'), &
       nctkarr_t('reduced_coordinates_plane_waves_dielectric_function', "i", &
       'number_of_reduced_dimensions, number_of_coefficients_dielectric_function, number_of_qpoints_dielectric_function'), &
       ! Abinit
       nctkarr_t('vcutgeo', "dp", 'number_of_reduced_dimensions'), &
       nctkarr_t('kind_cdata', "char", 'character_string_length'), &
       nctkarr_t('titles', "char", 'character_string_length, two') &
     ])
     NCF_CHECK(ncerr)

     ! FIXME problem with q-points, heads and wings?
     ! The order in P in the specs is wrong, q should be the last dimension here I use the "correct" version
     ! TODO: Remove. Use abifile_t
     varname = ncname_from_id(hscr%id)
     ncerr = nctk_def_arrays(ncid, &
       nctkarr_t(varname, "dp", &
&"complex, number_of_coefficients_dielectric_function, number_of_coefficients_dielectric_function,&
&number_of_spins, number_of_spins, number_of_frequencies_dielectric_function, number_of_qpoints_dielectric_function"))
     NCF_CHECK(ncerr)

     !write(std_out,*)"nqlwl",hscr%nqlwl
     NCF_CHECK(nctk_set_datamode(ncid))
     call c_f_pointer(c_loc(hscr%omega(1)), r2vals, shape=[2, size(hscr%omega)])
     NCF_CHECK(nf90_put_var(ncid, vid('frequencies_dielectric_function'), r2vals))
     NCF_CHECK(nf90_put_var(ncid, vid('qpoints_dielectric_function'), hscr%qibz))
     NCF_CHECK(nf90_put_var(ncid, vid('reduced_coordinates_plane_waves_dielectric_function'), hscr%gvec))

     NCF_CHECK(nf90_put_var(ncid, vid("titles"), hscr%titles))
     NCF_CHECK(nf90_put_var(ncid, vid("kind_cdata"), hscr%kind_cdata))
     NCF_CHECK(nf90_put_var(ncid, vid("vcutgeo"), hscr%vcutgeo))

     ncerr = nctk_defnwrite_ivars(ncid, [character(len=nctk_slen) :: &
       "id", "ikxc", "inclvkb", "headform", "fform", "gwcalctyp", &
       "nbands_used", "npwwfn_used", "spmeth", "test_type", "tordering", "awtr", "icutcoul", &
       "gwcomp", "gwgamma" &
      ],&
      [ hscr%id, hscr%ikxc, hscr%inclvkb, hscr%headform, hscr%fform, hscr%gwcalctyp, &
       hscr%nbnds_used, hscr%npwwfn_used, hscr%spmeth, hscr%test_type, hscr%tordering, hscr%awtr, hscr%icutcoul, &
       hscr%gwcomp, hscr%gwgamma &
     ])
     NCF_CHECK(ncerr)

     ncerr = nctk_defnwrite_dpvars(ncid, [character(len=nctk_slen) :: &
        "mbpt_sciss", "spsmear", "zcut", "gwencomp"], &
        [hscr%mbpt_sciss, hscr%spsmear, hscr%zcut, hscr%gwencomp &
     ])
     NCF_CHECK(ncerr)

     ! TODO
     ! Add q-points for heads and wings for q-->0.
     if (hscr%nqlwl>0) then
       head_shape = "complex, number_of_spins, number_of_spins, number_of_frequencies_dielectric_function"
       head_shape = trim(head_shape)//", number_of_qpoints_gamma_limit"

       wing_shape = "complex, number_of_coefficients_dielectric_function, number_of_spins, number_of_spins"
       wing_shape = trim(wing_shape)//", number_of_frequencies_dielectric_function, number_of_qpoints_gamma_limit"

       ncerr = nctk_def_arrays(ncid, [&
         nctkarr_t("dielectric_function_head", "dp", head_shape),&
         nctkarr_t("dielectric_function_upper_wing", "dp", wing_shape),&
         nctkarr_t("dielectric_function_lower_wing", "dp", wing_shape)], defmode=.True.)
       NCF_CHECK(ncerr)

       NCF_CHECK(nctk_set_datamode(ncid))
       NCF_CHECK(nf90_put_var(ncid, vid('qpoints_gamma_limit'), hscr%qlwl))
     end if
#endif
   else
     MSG_ERROR(sjoin('Unsupported iomode:',iomode2str(iomode)))
   end if

 else
   MSG_BUG(sjoin("Wrong value for rdwr:", itoa(rdwr)))
 end if ! read/write/echo

 DBG_EXIT("COLL")

 return

 ! Handle Fortran IO error
 10 continue
 MSG_ERROR(errmsg)

contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end subroutine hscr_io
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/hscr_print
!! NAME
!! hscr_print
!!
!! FUNCTION
!!  Prints info on the header of the SCR|SUSC file.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_io_screening,m_screen,m_screening,mrgscr,setup_bse
!!
!! CHILDREN
!!      hscr_copy,hscr_io,hscr_print,read_screening,write_screening,wrtout
!!
!! SOURCE

subroutine hscr_print(Hscr,header,unit,prtvol,mode_paral)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: prtvol,unit
 character(len=4),intent(in),optional :: mode_paral
 character(len=*),intent(in),optional :: header
 type(hscr_t),intent(in) :: hscr

!Local variables-------------------------------
!scalars
 integer :: iomega,iqibz,unt,verbose
 character(len=4) :: mode
 character(len=500) :: msg

! *************************************************************************

 unt=std_out; if (PRESENT(unit      )) unt    =unit
 verbose=0  ; if (PRESENT(prtvol    )) verbose=prtvol
 mode='COLL'; if (PRESENT(mode_paral)) mode   =mode_paral

 if (PRESENT(header)) then
   msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
   call wrtout(unt,msg,mode)
 end if

 write(msg,'(1x,a)')TRIM(hscr%titles(1))
 call wrtout(unt,msg,mode)
 write(msg,'(1x,a)')TRIM(hscr%titles(2))
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Identifier                ',hscr%ID
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Kxc kernel                ',hscr%ikxc
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Treatment of q-->0 limit  ',hscr%inclvkb
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' headform                  ',hscr%headform
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' fform                     ',hscr%fform
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' gwcalctyp                 ',hscr%gwcalctyp
 call wrtout(unt,msg,mode)
 write(msg,'(a,2i8)')' Number of components      ',hscr%nI,hscr%nJ
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Number of q-points        ',hscr%nqibz
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Number of q-directions    ',hscr%nqlwl
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Number of frequencies     ',hscr%nomega
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Number of bands used      ',hscr%nbnds_used
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Dimension of matrix       ',hscr%npwe
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Number of planewaves used ',hscr%npwwfn_used
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Spectral method           ',hscr%spmeth
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Test_type                 ',hscr%test_type
 call wrtout(unt,msg,mode)
 write(msg,'(a,i8)') ' Time-ordering             ',hscr%tordering
 call wrtout(unt,msg,mode)
 write(msg,'(a,es16.6)')' Scissor Energy             ',hscr%mbpt_sciss
 call wrtout(unt,msg,mode)
 write(msg,'(a,es16.6)')' Spectral smearing          ',hscr%spsmear
 call wrtout(unt,msg,mode)
 write(msg,'(a,es16.6)')' Complex Imaginary Shift    ',hscr%zcut
 call wrtout(unt,msg,mode)

 if (verbose==0) then
   call wrtout(unt,' The header contains additional records.',mode)
 else
   write(msg,'(2a)')ch10,' q-points [r.l.u.]:'
   call wrtout(unt,msg,mode)
   do iqibz=1,hscr%nqibz
     write(msg,'(i5,3f12.6)')iqibz,hscr%qibz(:,iqibz)
     call wrtout(unt,msg,mode)
   end do

   write(msg,'(2a)')ch10,' Frequencies used [eV]:'
   call wrtout(unt,msg,mode)
   do iomega=1,hscr%nomega
     write(msg,'(i3,2f7.2)')iomega,REAL(hscr%omega(iomega))*Ha_eV,AIMAG(hscr%omega(iomega))*Ha_eV
     call wrtout(unt,msg,mode)
   end do
 end if

! HSCR_NEW
! HSCR_NEW

 ! Echo the abinit header.
 !if (prtvol>0) call hdr_echo(hscr%hdr,fform,rdwr,unit=unt)

end subroutine hscr_print
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/hscr_new
!! NAME
!!  hscr_new
!!
!! FUNCTION
!!  Initialize the Hscr datatype and most of its content from the
!!  em1params_t data type Ep.
!!
!! INPUTS
!!  varname=Name of the netcdf variable (used to get fform and ID).
!!  ikxc=Integer flag definining the type of XC kernel (0 if None i.e RPA)
!!  test_type=Integer flag defining the type of probing charge (0 for None)
!!  tordering=The time-ordering of the Response function.
!!  gvec(3,Ep%npwe)=The G-vectors used.
!!  Ep<em1params_t>=Parameters defining the calculation of the screening.
!!  hdr_abinit<hdr_type>=The abinit header.
!!
!! OUTPUT
!!  Hscr<type(hscr_t)>=the header, initialized.
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!
!! SOURCE

type(hscr_t) function hscr_new(varname,dtset,ep,hdr_abinit,ikxc,test_type,tordering,titles,ngvec,gvec) result(hscr)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikxc,test_type,tordering,ngvec
 character(len=*),intent(in) :: varname
 type(dataset_type),intent(in) :: dtset
 type(em1params_t),intent(in) :: Ep
 type(hdr_type),intent(in) :: hdr_abinit
!arrays
 integer,intent(in) :: gvec(3,ngvec)
 character(len=80),intent(in) :: titles(2)

!Local variables-------------------------------
 integer :: id
 type(abifile_t) :: abifile
! *************************************************************************

 !@hscr_t
 ABI_CHECK(ngvec==Ep%npwe,'ngvec/=Ep%npwe')

 ! ID=Identifier used to define the type of Response function (e^-1, chi0)
 id = 0
 if (varname == "polarizability") id = 1
 if (varname == "inverse_dielectric_function") id = 4
 ABI_CHECK(id /= 0, sjoin("Invalid varname: ",varname))

 ! Get fform from abifile.
 abifile = abifile_from_varname(varname)
 if (abifile%fform == 0) then
    MSG_ERROR(sjoin("Cannot find any abifile object associated to varname:", varname))
 end if

 ! Copy the abinit header.
 call hdr_copy(hdr_abinit,Hscr%Hdr)

 ! Initialize quantities related to the screening file
 hscr%id         =id
 hscr%ikxc       =ikxc
 hscr%inclvkb    =Ep%inclvkb
 hscr%headform   =HSCR_LATEST_HEADFORM
 hscr%fform      =abifile%fform
 hscr%gwcalctyp  =Ep%gwcalctyp
 hscr%nI         =Ep%nI
 hscr%nJ         =Ep%nJ
 hscr%nqibz      =Ep%nqcalc  ! nqcalc == nqibz except if we split the calculation with nqptdm
 hscr%nqlwl      =Ep%nqlwl
 hscr%nomega     =Ep%nomega
 hscr%nbnds_used =Ep%nbnds
 hscr%npwe       =Ep%npwe
 hscr%npwwfn_used=Ep%npwwfn
 hscr%spmeth     =Ep%spmeth
 hscr%test_type  =test_type
 hscr%tordering  =tordering
 hscr%mbpt_sciss =Ep%mbpt_sciss
 hscr%spsmear    =Ep%spsmear
 hscr%zcut       =Ep%zcut

 hscr%titles(:)=titles(:)

 call hscr_malloc(hscr, hscr%npwe, hscr%nqibz, hscr%nomega, hscr%nqlwl)
 hscr%gvec(:,:) = gvec(1:3,1:Ep%npwe)
 hscr%qibz(:,:) = Ep%qcalc
 hscr%qlwl(:,:) = Ep%qlwl
 hscr%omega(:) = Ep%omega

! HSCR_NEW
 hscr%awtr = dtset%awtr
 hscr%icutcoul = dtset%icutcoul
 hscr%vcutgeo = dtset%vcutgeo
 hscr%gwcomp = dtset%gwcomp
 hscr%gwgamma = dtset%gwgamma
 hscr%gwencomp = dtset%gwencomp
 hscr%kind_cdata = "dpc" ! For the time being, data is always written in double-precision
! HSCR_NEW

end function hscr_new
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/hscr_bcast
!! NAME
!! hscr_bcast
!!
!! FUNCTION
!! This subroutine transmit the header structured datatype initialized
!! on one processor (or a group of processor), to the other processors.
!! It also allocates the needed part of the header.
!!
!! INPUTS
!!  master=ID of the master node.
!!  my_rank=ID of the node that receives the data.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  (no output)
!!
!! SIDE EFFECTS
!!  Hscr<type(hscr_t)>=the SCR header. For the master, it is already
!!   initialized entirely, while for the other procs, everything has
!!   to be transmitted.
!!
!! NOTES
!! This routine is called only in the case of MPI version of the code.
!!
!! PARENTS
!!      m_io_screening,setup_bse
!!
!! CHILDREN
!!      hscr_copy,hscr_io,hscr_print,read_screening,write_screening,wrtout
!!
!! SOURCE

subroutine hscr_bcast(hscr,master,my_rank,comm)

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: master,my_rank,comm
 type(hscr_t),intent(inout) :: hscr

!Local variables-------------------------------
 integer :: ierr

! *************************************************************************

 DBG_ENTER("COLL")
 if (xmpi_comm_size(comm) == 1) return ! Nothing to do

 !@hscr_t
! integer
 call xmpi_bcast(hscr%id,         master,comm,ierr)
 call xmpi_bcast(hscr%ikxc,       master,comm,ierr)
 call xmpi_bcast(hscr%inclvkb,    master,comm,ierr)
 call xmpi_bcast(hscr%headform,   master,comm,ierr)
 call xmpi_bcast(hscr%fform,      master,comm,ierr)
 call xmpi_bcast(hscr%gwcalctyp,  master,comm,ierr)
 call xmpi_bcast(hscr%nI,         master,comm,ierr)
 call xmpi_bcast(hscr%nJ,         master,comm,ierr)
 call xmpi_bcast(hscr%nqibz,      master,comm,ierr)
 call xmpi_bcast(hscr%nqlwl,      master,comm,ierr)
 call xmpi_bcast(hscr%nomega,     master,comm,ierr)
 call xmpi_bcast(hscr%nbnds_used, master,comm,ierr)
 call xmpi_bcast(hscr%npwe,       master,comm,ierr)
 call xmpi_bcast(hscr%npwwfn_used,master,comm,ierr)
 call xmpi_bcast(hscr%spmeth,     master,comm,ierr)
 call xmpi_bcast(hscr%test_type,  master,comm,ierr)
 call xmpi_bcast(hscr%tordering,  master,comm,ierr)

 ! Real
 call xmpi_bcast(hscr%mbpt_sciss, master,comm,ierr)
 call xmpi_bcast(hscr%spsmear,    master,comm,ierr)
 call xmpi_bcast(hscr%zcut,       master,comm,ierr)

 ! arrays
 call xmpi_bcast(hscr%titles, master,comm,ierr)

 if (my_rank /= master) then
   call hscr_malloc(hscr, hscr%npwe, hscr%nqibz, hscr%nomega, hscr%nqlwl)
 end if

 call xmpi_bcast(hscr%gvec, master,comm,ierr)
 call xmpi_bcast(hscr%qibz, master,comm,ierr)
 call xmpi_bcast(hscr%qlwl, master,comm,ierr)
 call xmpi_bcast(hscr%omega,master,comm,ierr)

 ! Communicate the Abinit header.
 call hdr_bcast(hscr%Hdr,master,my_rank,comm)

! HSCR_NEW
 call xmpi_bcast(hscr%awtr, master, comm, ierr)
 call xmpi_bcast(hscr%icutcoul, master, comm, ierr)
 call xmpi_bcast(hscr%vcutgeo, master, comm, ierr)
 call xmpi_bcast(hscr%gwcomp, master, comm, ierr)
 call xmpi_bcast(hscr%gwgamma, master, comm, ierr)
 call xmpi_bcast(hscr%gwencomp, master, comm, ierr)
 call xmpi_bcast(hscr%kind_cdata, master, comm, ierr)
! HSCR_NEW

 DBG_EXIT("COLL")

end subroutine hscr_bcast
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/hscr_malloc
!! NAME
!! hscr_malloc
!!
!! FUNCTION
!! Allocate the components of the header structured datatype except for hscr%hdr
!!
!! PARENTS
!!      m_io_screening
!!
!! CHILDREN
!!      hscr_copy,hscr_io,hscr_print,read_screening,write_screening,wrtout
!!
!! SOURCE

subroutine hscr_malloc(hscr, npwe, nqibz, nomega, nqlwl)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwe, nqibz, nomega, nqlwl
 type(hscr_t),intent(inout) :: Hscr

! *************************************************************************

 !@hscr_t
 ABI_MALLOC(hscr%gvec, (3, npwe))
 ABI_MALLOC(hscr%qibz, (3, nqibz))
 ABI_MALLOC(hscr%qlwl, (3, nqlwl))
 ABI_MALLOC(hscr%omega, (nomega))

end subroutine hscr_malloc
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/hscr_free
!! NAME
!! hscr_free
!!
!! FUNCTION
!! Deallocate the components of the header structured datatype
!!
!! INPUTS
!! hdr <type(hdr_type)>=the header
!!
!! OUTPUT
!!  (only deallocate)
!!
!! PARENTS
!!      m_io_screening,m_screen,m_screening,mrgscr,screening,setup_bse
!!
!! CHILDREN
!!      hscr_copy,hscr_io,hscr_print,read_screening,write_screening,wrtout
!!
!! SOURCE

subroutine hscr_free(hscr)

 implicit none

!Arguments ------------------------------------
!scalars
 type(hscr_t),intent(inout) :: hscr

! *************************************************************************

 !@hscr_t
 DBG_ENTER("COLL")

 if (allocated(hscr%gvec)) then
   ABI_FREE(hscr%gvec)
 end if
 if (allocated(hscr%qibz)) then
   ABI_FREE(hscr%qibz)
 end if
 if (allocated(hscr%qlwl)) then
   ABI_FREE(hscr%qlwl)
 end if
 if (allocated(hscr%omega)) then
   ABI_FREE(hscr%omega)
 end if

 call hdr_free(hscr%Hdr)

 DBG_EXIT("COLL")

end subroutine hscr_free
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/hscr_copy
!! NAME
!! hscr_copy
!!
!! FUNCTION
!! Deep copy of the header of the _SCR or _SUSC file.
!!
!! INPUTS
!!
!! PARENTS
!!      m_io_screening,m_screening
!!
!! CHILDREN
!!      hscr_copy,hscr_io,hscr_print,read_screening,write_screening,wrtout
!!
!! SOURCE

subroutine hscr_copy(Hscr_in,Hscr_cp)

 implicit none

!Arguments ------------------------------------
!scalars
 type(hscr_t),intent(in) :: Hscr_in
 type(hscr_t),intent(inout) :: Hscr_cp

!Local variables-------------------------------
!scalars
 !character(len=500) :: msg
! *************************************************************************

 !@hscr_t
 ! Integer values.
 Hscr_cp%id          = Hscr_in%id
 Hscr_cp%ikxc        = Hscr_in%ikxc
 Hscr_cp%inclvkb     = Hscr_in%inclvkb
 Hscr_cp%headform    = Hscr_in%headform
 Hscr_cp%fform       = Hscr_in%fform
 Hscr_cp%gwcalctyp   = Hscr_in%gwcalctyp
 Hscr_cp%nI          = Hscr_in%nI
 Hscr_cp%nJ          = Hscr_in%nJ
 Hscr_cp%nqibz       = Hscr_in%nqibz
 Hscr_cp%nqlwl       = Hscr_in%nqlwl
 Hscr_cp%nomega      = Hscr_in%nomega
 Hscr_cp%nbnds_used  = Hscr_in%nbnds_used
 Hscr_cp%npwe        = Hscr_in%npwe
 Hscr_cp%npwwfn_used = Hscr_in%npwwfn_used
 Hscr_cp%spmeth      = Hscr_in%spmeth
 Hscr_cp%test_type   = Hscr_in%test_type
 Hscr_cp%tordering   = Hscr_in%tordering

 ! Real variables
 Hscr_cp%mbpt_sciss = Hscr_in%mbpt_sciss
 Hscr_cp%spsmear  = Hscr_in%spsmear
 Hscr_cp%zcut     = Hscr_in%zcut

 ! Copy the abinit Header
 call hdr_copy(Hscr_in%Hdr,Hscr_cp%Hdr)

 Hscr_cp%titles(:) = Hscr_in%titles(:)

 ! Copy allocatable arrays.
 call alloc_copy(Hscr_in%gvec , Hscr_cp%gvec)
 call alloc_copy(Hscr_in%qibz , Hscr_cp%qibz)
 call alloc_copy(Hscr_in%qlwl , Hscr_cp%qlwl)
 call alloc_copy(Hscr_in%omega, Hscr_cp%omega)

! HSCR_NEW
 hscr_cp%awtr      =  hscr_in%awtr
 hscr_cp%icutcoul  =  hscr_in%icutcoul
 hscr_cp%vcutgeo   =  hscr_in%vcutgeo
 hscr_cp%gwcomp    =  hscr_in%gwcomp
 hscr_cp%gwgamma   =  hscr_in%gwgamma
 hscr_cp%gwencomp  =  hscr_in%gwencomp
 hscr_cp%kind_cdata  =  hscr_in%kind_cdata
! HSCR_NEW

end subroutine hscr_copy
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/hscr_merge
!! NAME
!! hscr_merge
!!
!! FUNCTION
!! This subroutine merges diffrent header structured variable (hscr_t)
!!
!! INPUTS
!!  Hscr_in(:) <hscr_t)>=List of headers to be merged.
!!
!! OUTPUT
!!  Hscr_out<hscr_t>=The output merged header.
!!
!! PARENTS
!!      m_io_screening
!!
!! CHILDREN
!!      hscr_copy,hscr_io,hscr_print,read_screening,write_screening,wrtout
!!
!! SOURCE

subroutine hscr_merge(Hscr_in,Hscr_out)

 implicit none

!Arguments ------------------------------------
!scalars
 type(hscr_t),intent(in) :: Hscr_in(:)
 type(hscr_t),intent(out) :: Hscr_out

!Local variables-------------------------------
!scalars
 integer :: nhds,restart,restartpaw,ihd,ii,nqtot,nqneq
 logical :: isok
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: qset(:,:)
! *************************************************************************

 !@hscr_t
 nhds=SIZE(Hscr_in)

 ! Initial copy of the header ===
 ! If multiple headers, select the header containing q-->0 so that we copy also heads and wings
 ii = imax_loc(Hscr_in(:)%nqlwl)
 call hscr_copy(Hscr_in(ii),Hscr_out)
 if (nhds==1) return

 ! Check consistency of the abinit Headers.
 ! FFT grid might be q-point dependent so we stop only when restart==0
 isok=.TRUE.
 do ihd=2,nhds
   call hdr_check(Hscr_in(1)%fform,Hscr_in(ihd)%fform,Hscr_in(1)%Hdr,Hscr_in(ihd)%Hdr,'COLL',restart,restartpaw)
   if (restart==0) then
     isok=.FALSE.
     write(msg,'(a,i0,a)')' Abinit header no.',ihd,' is not consistent with the first header '
     MSG_WARNING(msg)
   end if
 end do
 if (.not.isok) then
   MSG_ERROR('Cannot continue, Check headers')
 end if

 ! Now check variables related to polarizability|epsilon^{-1}.
 ! 1) Tests quantities that must be equal
 ii = assert_eq(Hscr_in(:)%ID,       'Headers have different Identifiers')
 ii = assert_eq(Hscr_in(:)%ikxc,     'Headers have different ikxc'       )
 ii = assert_eq(Hscr_in(:)%headform, 'Headers have different headform'   )
 ii = assert_eq(Hscr_in(:)%fform,    'Headers have different fform'      )
 ii = assert_eq(Hscr_in(:)%gwcalctyp,'Headers have different gwcalctyp'  )
 ii = assert_eq(Hscr_in(:)%nI,       'Headers have different nI'         )
 ii = assert_eq(Hscr_in(:)%nJ,       'Headers have different nJ'         )
 ii = assert_eq(Hscr_in(:)%nomega,   'Headers have different nomega'     )
 ii = assert_eq(Hscr_in(:)%test_type,'Headers have different test_type'  )
 ii = assert_eq(Hscr_in(:)%tordering,'Headers have different tordering'  )

 ! This is not mandatory but makes life easier!
 ii = assert_eq(Hscr_in(:)%npwe,'Headers have different number of G-vectors'  )

 do ihd=2,nhds
   if (ANY(ABS(Hscr_in(ihd)%omega-Hscr_in(1)%omega)>tol6)) then
     write(msg,'(a,i0,a)')' Frequencies in the first and the ',ihd,'-th header differ'
     MSG_ERROR(msg)
   end if
   if (ANY(Hscr_in(ihd)%gvec(:,:)-Hscr_in(1)%gvec(:,:)/=0)) then
     write(msg,'(a,i0,a)')' Incompatible G-vector list found in the ',ihd,'-th header'
     MSG_ERROR(msg)
   end if
   if (hscr_in(ihd)%kind_cdata /= hscr_in(1)%kind_cdata) then
     write(msg,'(3a,i0,2a)')' Files contain data with different precisions.',ch10,&
     "In particular the ",ihd,'-th header has precision:',trim(hscr_in(ihd)%kind_cdata)
     MSG_ERROR(msg)
   end if
 end do !ihd

 ! If error is not fatal, just warn ===
 if (ANY(Hscr_in(:)%npwwfn_used/=Hscr_in(1)%npwwfn_used)) then
   MSG_COMMENT('Files have been produced with a different number of planewaves for the wavefunctions.')
 end if
 if (ANY(Hscr_in(:)%nbnds_used/=Hscr_in(1)%nbnds_used)) then
   MSG_COMMENT('Files have been produced with a different number of bands.')
 end if
 if (ANY(Hscr_in(:)%spmeth/=Hscr_in(1)%spmeth)) then
   MSG_COMMENT('Files have been produced with different algorithms.')
 end if
 if (ANY(ABS(Hscr_in(:)%mbpt_sciss-Hscr_in(1)%mbpt_sciss)>tol6)) then
   MSG_COMMENT('Files have benn produced with different values of mbpt_sciss.')
 end if
 if (ANY(ABS(Hscr_in(:)%spsmear-Hscr_in(1)%spsmear)>tol6)) then
   MSG_COMMENT('Files have been produced with different values of spsmear.')
 end if
 if (ANY(ABS(Hscr_in(:)%zcut-Hscr_in(1)%zcut)>tol6)) then
   MSG_COMMENT('Files have been produced with different values of zcut.')
 end if

 ! Now merge the list of q-points.
 ! Take the union of the q-points, remove possible duplicated
 ! are change the parameters in hscr_out that depends on q-points.
 nqtot=SUM(Hscr_in(:)%nqibz)
 ABI_MALLOC(qset,(3,nqtot))

 ii=0
 do ihd=1,nhds
   qset(:,ii+1:ii+Hscr_in(ihd)%nqibz)=Hscr_in(ihd)%qibz(:,:)
   ii=ii+Hscr_in(ihd)%nqibz
 end do

 call remove_copies(nqtot,qset,nqneq,isequalk)

 if (nqneq /= nqtot) then
   write(msg,'(3a,2(i0,a))')&
    'COMMENT: Headers contain duplicated q-points ',ch10,&
    'Found ',nqneq,' distinct q-points among the total ',nqtot,' points reported in the headers. '
   call wrtout(std_out, msg)
 end if

 Hscr_out%nqibz = nqneq
 ABI_FREE(Hscr_out%qibz)
 ABI_MALLOC(Hscr_out%qibz,(3,nqneq))
 Hscr_out%qibz(:,:)=qset(:,1:nqneq)
 ABI_FREE(qset)

end subroutine hscr_merge
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/write_screening
!! NAME
!! write_screening
!!
!! FUNCTION
!! For a single q-point, write either \tilde epsilon^{-1} on the _SCR file
!! or chi0 on the _SUSC file. The file is supposed to have been open in the calling routine.
!!
!! INPUTS
!!  varname=The name of the array to write (used if etsf-io format).
!!  unt=The unit number of the file to be written (supposed to be already open)
!!  iomode=Integer flag defining the format of the output file. Available options:
!!    IO_MODE_FORTRAN--> Plain Fortran file
!!    IO_MODE_ETSF--> ETSF format
!!  npwe=Number of plane waves in epsm1.
!!  nomega=Number of frequencies
!!  iqibz=Index of the q-points in the IBZ.
!!  epsm1(npwe,npwe,nomega)=The matrix to be written, for different frequencies, and a single q-point.
!!
!! NOTES
!!  On some architecture, the code crashes when trying to write or read a record containing the
!!  entire (G1,G2) matrix thus we use smaller records containing the columns of the two-point function.
!!
!! OUTPUT
!!  (only writing on file)
!!
!! PARENTS
!!      m_io_screening,m_screening,screening
!!
!! CHILDREN
!!      hscr_copy,hscr_io,hscr_print,read_screening,write_screening,wrtout
!!
!! SOURCE

subroutine write_screening(varname,unt,iomode,npwe,nomega,iqibz,epsm1)

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: varname
 integer,intent(in) :: nomega,npwe,iqibz,unt,iomode
!arrays
 complex(gwpc),target,intent(in) :: epsm1(npwe,npwe,nomega)

!Local variables-------------------------------
!scalars
 integer :: ipwe,iomega,spins(2),s1,s2
 character(len=500) :: errmsg
!arrays
 complex(dpc),allocatable :: epsm1d(:,:)
#ifdef HAVE_NETCDF
 integer :: varid,ncerr
#endif
#ifdef HAVE_GW_DPC
 real(dp), ABI_CONTIGUOUS pointer :: real_epsm1(:,:,:,:,:,:,:)
#else
 real(sp), ABI_CONTIGUOUS pointer :: real_epsm1(:,:,:,:,:,:,:)
#endif

! *************************************************************************

 DBG_ENTER("COLL")

 select case (iomode)
 case (IO_MODE_FORTRAN, IO_MODE_MPI)
   ! Write a record for each omega, Always use double precision.
   ABI_MALLOC(epsm1d,(npwe,1))

   do iomega=1,nomega
     do ipwe=1,npwe
       epsm1d(:,1) = epsm1(:,ipwe,iomega) !spc ==> dpc
       write(unt, err=10, iomsg=errmsg)epsm1d(1:npwe,1)
     end do
   end do
   ABI_FREE(epsm1d)

#ifdef HAVE_NETCDF
 case (IO_MODE_ETSF)
   ! netcdf does not support complex datatypes. Here I use some C-magic to  associate the memory
   ! to a Fortran real pointer with the correct type and shape. Note that the data on file is always in double precision.
   varid = nctk_idname(unt, varname)
   call c_f_pointer(c_loc(epsm1(1,1,1)), real_epsm1, [2,npwe,npwe,1,1,nomega,1])
   ! [cplex, npwe, npwe, nspin, nspin, nomega, nqpt]
   spins = 1; s1 = spins(1); s2 = spins(2)
   ncerr = nf90_put_var(unt, varid, real_epsm1, start=[1,1,1,s1,s2,1,iqibz], count=[2,npwe,npwe,1,1,nomega,1])
   NCF_CHECK_MSG(ncerr, sjoin("putting var:", varname))
#endif

 case default
   MSG_ERROR(sjoin("Wrong iomode:", iomode2str(iomode)))
 end select

 DBG_EXIT("COLL")

 return

 ! Handle IO error
 10 continue
 MSG_ERROR(errmsg)

end subroutine write_screening
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/read_screening
!! NAME
!! read_screening
!!
!! FUNCTION
!! Read either a screening (\tilde epsilon^{-1}) file in the SCR format or
!! the irreducible polarizability (chi0) in the SUSC format.
!!
!! INPUTS
!!  varname=Name of the array to read. Used for ETSF-IO files.
!!  iomode=Integer flag defining the format of the output file. Available options:
!!    IO_MODE_FORTRAN--> Plain Fortran file
!!    IO_MODE_ETSF--> ETSF format
!!  iqiA[optional]=Used if only a particular q-point is required. In this case iqiA define the index
!!   of the required q-point in the array qibz(3,Hscr%nqibz)
!!  nqibzA=number of asked q-points (used to dimension the output arrays).
!!   Equal to Hscr%nqibz if the full matrix is required
!!  comm=MPI communicator.
!!  npweA=number of asked planewaves
!!  nomegaA=number of asked frequencies
!!
!! OUTPUT
!!  epsm1(npweA,npweA,nomegaA,nqibzA) = \tilde\epsilon^{-1}(Ng,Ng,Nw,Nq)
!!
!! NOTES
!!  * If the epsilon matrix read is bigger than npweA x npweA, it will be truncated;
!!    if it is smaller, an error will occur
!!  * If the number of frequencies asked for is smaller than that reported in the file, the matrix
!!    will be truncated. If nomegaA > Hscr%nomega an error will occur
!!
!! PARENTS
!!      calc_ucrpa,m_io_screening,m_screen,m_screening,mrgscr
!!
!! CHILDREN
!!      hscr_copy,hscr_io,hscr_print,read_screening,write_screening,wrtout
!!
!! SOURCE

subroutine read_screening(varname,fname,npweA,nqibzA,nomegaA,epsm1,iomode,comm,&
& iqiA) ! Optional

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iomode,nomegaA,npweA,nqibzA,comm
 integer,optional,intent(in) :: iqiA
 character(len=*),intent(in) :: varname,fname
!arrays
 complex(gwpc),target,intent(inout) :: epsm1(npweA,npweA,nomegaA,nqibzA)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: ipwe,fform,iomega,iqibz,unt,rdwr,my_rank,nprocs,my_iomode
#ifdef HAVE_NETCDF
 integer :: varid,ncerr
#endif
#ifdef HAVE_MPI_IO
 integer :: test_fform,mpi_err,ierr,sc_mode
 integer :: bsize_frm,mpi_type_frm
 integer :: mpi_fh,buf_dim !,mat_ggw,mat_ggwq
 integer(XMPI_OFFSET_KIND) :: offset,displ_wq !,my_offpad
 !complex(dpc) :: ctmp
#endif
 character(len=500) :: msg,errmsg
 logical :: read_qslice
 type(hscr_t) :: Hscr
!arrays
#ifdef HAVE_MPI_IO
 integer(MPI_OFFSET_KIND),allocatable :: offset_wq(:,:)
#endif
 complex(dpc),allocatable :: bufdc2d(:,:),bufdc3d(:,:,:)
 ! pointers passed to netcdf4 routines (complex datatypes are not supported).
#ifdef HAVE_GW_DPC
 real(dp), ABI_CONTIGUOUS pointer :: real_epsm1(:,:,:,:,:,:,:)
#else
 real(sp), ABI_CONTIGUOUS pointer :: real_epsm1(:,:,:,:,:,:,:)
#endif
 integer :: spins(2),s1,s2

! *************************************************************************

 DBG_ENTER("COLL")

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 my_iomode = iomode
 if (endswith(fname, ".nc")) my_iomode = IO_MODE_ETSF
 !my_iomode = IO_MODE_MPI
 !if (my_iomode  == IO_MODE_MPI) my_iomode = IO_MODE_FORTRAN

 rdwr=1
 select case (my_iomode)
 case (IO_MODE_MPI)
#ifdef HAVE_MPI_IO
   bsize_frm    = xmpio_bsize_frm    ! bsize_frm= Byte length of the Fortran record marker.
   mpi_type_frm = xmpio_mpi_type_frm ! MPI type of the record marker.
   sc_mode = xmpio_collective

   ! Master reads the header via Fortran IO then bcast the data.
   call hscr_from_file(hscr, fname, fform, comm)

   ! Open the file with MPI-IO
   call MPI_FILE_OPEN(comm, fname, MPI_MODE_RDONLY, xmpio_info ,mpi_fh, mpi_err)
   ABI_CHECK_MPI(mpi_err, sjoin("MPI_FILE_OPEN:", fname))

   ! Retrieve the offset of the section immediately below the header.
   call hscr_mpio_skip(mpi_fh,test_fform,offset)
   ABI_CHECK(test_fform == fform, "mismatch in fform!")

   ! Offsets of the Fortran markers corresponding to the (w,q) slices.
   ABI_MALLOC(offset_wq,(HScr%nomega,HScr%nqibz))
   displ_wq = offset
   do iqibz=1,Hscr%nqibz
     do iomega=1,Hscr%nomega
       ABI_CHECK(displ_wq > 0, "displ_wq < 0, your SCR|SUSC file is too big for MPI-IO!")
       offset_wq(iomega,iqibz) = displ_wq
       displ_wq = displ_wq + Hscr%npwe**2 * xmpi_bsize_dpc + Hscr%npwe * 2 * bsize_frm
     end do
   end do
#else
   MSG_ERROR("MPI-IO support not enabled at configure-time")
#endif

 case (IO_MODE_FORTRAN)
   ! Plain Fortran IO, all nodes read.
   if (open_file(fname,msg,newunit=unt,form="unformatted",status="old",action="read") /= 0) then
     MSG_ERROR(msg)
   end if
   call hscr_io(hscr,fform,rdwr,unt,comm,master,my_iomode)

 case (IO_MODE_ETSF)
#ifdef HAVE_NETCDF
   NCF_CHECK(nctk_open_read(unt, fname, xmpi_comm_self))
   call hscr_io(hscr,fform,rdwr,unt,comm,master,my_iomode)
#endif

 case default
   MSG_ERROR(sjoin("Wrong iomode:", iomode2str(my_iomode)))
 end select

 ! Slice or full array?
 read_qslice = .False.
 if (PRESENT(iqiA)) then
   read_qslice = .True.
   !call wrtout(std_out, sjoin('. Reading q-slice for iq = ',itoa(iqiA),' from: ', fname))
   if (iqiA <= 0 .or. iqiA > Hscr%nqibz) then
     MSG_BUG('iqiA out of range')
   end if
 end if

 ! Do some check
 if (Hscr%npwe>npweA) then
   write(msg,'(a,i0,2a,i0)')&
&    'Total number of G-vectors reported on file = ',Hscr%npwe,ch10,&
&    'Reading a smaller matrix of dimension      = ',npweA
   MSG_COMMENT(msg)
 end if

 if (npweA>Hscr%npwe) then
   write(msg,'(2(a,i0))')' Dimension of matrix = ',Hscr%npwe," requiring a too big matrix = ",npweA
   MSG_ERROR(msg)
 end if

 ABI_CHECK(nqibzA  <= Hscr%nqibz, 'Requiring too many q-points')
 ABI_CHECK(nomegaA <= Hscr%nomega,'Requiring too many frequencies')

 select case (my_iomode)
 case (IO_MODE_MPI)
#ifdef HAVE_MPI_IO
   if (read_qslice) then
      call wrtout(std_out, "calling mpiotk to read_qslice", "COLL")
      buf_dim = (npweA)**2 * nomegaA
      offset = offset_wq(1,iqiA)
      sc_mode = xmpio_collective

#ifdef HAVE_GW_DPC
     ! Read in-place.
     call mpiotk_read_fsuba_dpc3D(mpi_fh,offset, [HScr%npwe,HScr%npwe,HScr%nomega], [npweA,npweA,nomegaA], [1,1,1],&
        buf_dim,epsm1,xmpio_chunk_bsize,sc_mode,comm,ierr)
     ABI_CHECK(ierr==0,"Fortran matrix too big")
#else
     ! Have to allocate workspace for dpc data.
     ! FIXME: Change the file format of the SCR and SUC file so that
     ! they are written in single precision if not HAVE_GW_DPC
     ABI_MALLOC_OR_DIE(bufdc3d,(npweA,npweA,nomegaA), ierr)

     call mpiotk_read_fsuba_dpc3D(mpi_fh,offset, [HScr%npwe,HScr%npwe,HScr%nomega], [npweA,npweA,nomegaA], [1,1,1],&
        buf_dim,bufdc3d,xmpio_chunk_bsize,sc_mode,comm,ierr)
     ABI_CHECK(ierr==0,"Fortran matrix too big")

     epsm1(:,:,:,1) = bufdc3d
     ABI_FREE(bufdc3d)
#endif

   else
     ! Full matrix (G,G',w,q) is needed.
     call wrtout(std_out, "calling mpiotk: Full matrix (G,G',w,q) is needed.")

#ifdef HAVE_GW_DPC
     ! Can read all data at once.
     buf_dim = (npweA)**2 * nomegaA * HScr%nqibz
     offset = offset_wq(1,1)
     sc_mode = xmpio_collective

     call mpiotk_read_fsuba_dpc4D(mpi_fh,offset,&
       [HScr%npwe,HScr%npwe,HScr%nomega,HScr%nqibz], [npweA,npweA,nomegaA,HScr%nqibz], [1,1,1,1],&
       buf_dim,epsm1,xmpio_chunk_bsize,sc_mode,comm,ierr)
     ABI_CHECK(ierr==0,"Fortran record too big")
#else
     ! Have to allocate workspace for dpc data.
     ABI_MALLOC_OR_DIE(bufdc3d,(npweA,npweA,nomegaA), ierr)
     sc_mode = xmpio_collective

     do iqibz=1,Hscr%nqibz
       offset = offset_wq(1,iqibz)
       buf_dim = (2*npweA)**2 * nomegaA

       call mpiotk_read_fsuba_dpc3D(mpi_fh,offset, &
        [HScr%npwe,HScr%npwe,HScr%nomega], [npweA,npweA,nomegaA], [1,1,1],&
         buf_dim,bufdc3d,xmpio_chunk_bsize,sc_mode,comm,ierr)
       ABI_CHECK(ierr==0,"Fortran matrix too big")

       epsm1(:,:,:,iqibz) = bufdc3d
     end do

     ABI_FREE(bufdc3d)
#endif
   end if

   call MPI_FILE_CLOSE(mpi_fh,mpi_err)
   ABI_FREE(offset_wq)
#endif

 case (IO_MODE_FORTRAN)
   ! Read epsilon^-1 with Fortran IO
   ! Allocate a single column to save memory.
   ! TODO re-merge the two cases.
   ABI_MALLOC(bufdc2d,(Hscr%npwe,1))

   ! Two coding for different case just to keep it readable.
   select case (read_qslice)
   case (.True.)
     ! Read only a slice of the full array (useful if the entire array is huge).
     !if (dim_wings==1) STOP 'not implemented'
     !TODO this has to be done in a cleaner way.
     qread_loop: &
&    do iqibz=1,Hscr%nqibz
       if (iqibz==iqiA) then
         do iomega=1,nomegaA
           do ipwe=1,Hscr%npwe
             read(unt, err=10, iomsg=errmsg) bufdc2d(1:Hscr%npwe,1)
             if (ipwe<=npweA) epsm1(1:npweA,ipwe,iomega,1)=bufdc2d(1:npweA,1)
           end do
         end do
         EXIT qread_loop ! Got data. Do not need to read file till the end.
       else
         ! Skip other q-points i.e bufdc2d(1:Hscr%npwe,1:Hscr%npwe)
         do iomega=1,Hscr%nomega
           do ipwe=1,Hscr%npwe
            read(unt, err=10, iomsg=errmsg)
           end do
         end do
       end if ! iqibz==iqiA
     end do qread_loop ! iqibz

   case (.False.)
     ! Read the entire array.
     do iqibz=1,Hscr%nqibz
       do iomega=1,nomegaA
         do ipwe=1,Hscr%npwe
          read(unt, err=10, iomsg=errmsg) bufdc2d(1:Hscr%npwe,1)
          if (ipwe<=npweA) epsm1(1:npweA,ipwe,iomega,iqibz)=bufdc2d(1:npweA,1)
         end do
       end do
       ! Skip other frequencies
       do iomega=nomegaA+1,Hscr%nomega
         do ipwe=1,Hscr%npwe
           read(unt, err=10, iomsg=errmsg)
         end do
       end do
     end do !iqibz
   end select

   close(unt)

 case (IO_MODE_ETSF)
#ifdef HAVE_NETCDF
   ! netcdf does not support complex datatypes. Here I use some C-magic to  associate the memory
   ! to a Fortran real pointer with the correct type and shape. Note that the data on file is always in double precision.
   ! nf90_get_var will automatically convert from double to single if the GW code is in single precision mode.
   ! This is the reason why I'm using CPP option in the declaration of real_epsm1.

   ! FIXME: Need to know the type to read
   !write(std_out,*)"in read_screening"
   varid = nctk_idname(unt, varname)

   ! [cplex, npwe, npwe, nspin, nspin, nomega, nqpt]
   call c_f_pointer(c_loc(epsm1(1,1,1,1)), real_epsm1, [2,npweA,npweA,1,1,nomegaA,nqibzA])
   spins = 1; s1 = spins(1); s2 = spins(2)
   if (read_qslice) then
     ncerr = nf90_get_var(unt, varid, real_epsm1, start=[1,1,1,s1,s2,1,iqia], count=[2,npweA,npweA,1,1,nomegaA,1])
   else
     ncerr = nf90_get_var(unt, varid, real_epsm1, start=[1,1,1,s1,s2,1,1], count=[2,npweA,npweA,1,1,nomegaA,nqibzA])
     !do iqibz=1,nqibzA
     !  !write(*,*)"epsm1: in read ",iqibz,epsm1(1:3,1,1,iqibz)
     !end do
   end if
   NCF_CHECK_MSG(ncerr, sjoin("getting var:", varname))
   NCF_CHECK(nf90_close(unt))
   !write(std_out,*)"read_screening done"
#endif

 case default
   MSG_ERROR(sjoin("Wrong iomode:", iomode2str(my_iomode)))
 end select

 ! Free memory
 if (allocated(bufdc2d)) then
   ABI_FREE(bufdc2d)
 end if
 if (allocated(bufdc3d)) then
   ABI_FREE(bufdc3d)
 end if

 call hscr_free(Hscr)

 DBG_EXIT("COLL")

 return

 ! Handle Fortran IO error.
10 continue
 MSG_ERROR(errmsg)

end subroutine read_screening
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/hscr_mpio_skip
!! NAME
!!  hscr_mpio_skip
!!
!! FUNCTION
!!   Skip the header of the (SCR|SUSC) file in MPI-IO mode. This routine uses local MPI-IO calls hence
!!   it can be safely called by master node only. Note however that in this case the
!!   offset has to be communicated to the other nodes.
!!
!! INPUTS
!!  mpio_fh=MPI-IO file handler
!!  fmarker_bsize   = Byte length of Fortran record marker.
!!  fmarker_mpi_type= MPI type of the Fortran record marker
!!
!! OUTPUT
!!  fform=kind of the array in the file
!!  offset=The offset of the Fortran record located immediately below the Abinit header.
!!
!! PARENTS
!!      m_io_screening
!!
!! CHILDREN
!!      hscr_copy,hscr_io,hscr_print,read_screening,write_screening,wrtout
!!
!! SOURCE

subroutine hscr_mpio_skip(mpio_fh,fform,offset)

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: mpio_fh
 integer,intent(out) :: fform
 integer(kind=XMPI_OFFSET_KIND),intent(out) :: offset

!Local variables-------------------------------
!scalars
 integer :: bsize_frm,mpi_type_frm
#ifdef HAVE_MPI_IO
 integer :: ierr,isk,nqlwl
 !character(len=500) :: msg
!arrays
 integer(kind=MPI_OFFSET_KIND) :: fmarker,positloc
 integer :: statux(MPI_STATUS_SIZE)
#endif

! *************************************************************************

 offset = 0
 bsize_frm    = xmpio_bsize_frm    ! Byte size of the Fortran record marker.
 mpi_type_frm = xmpio_mpi_type_frm ! MPI type of the record marker.

 call hdr_mpio_skip(mpio_fh,fform,offset)

 call wrtout(std_out, sjoin("in hdr_mpio_skip with fform = ",itoa(fform)))

#ifdef HAVE_MPI_IO
 select case (fform)
 case (1003, 1004)
   ! Ship the titles
   call xmpio_read_frm(mpio_fh,offset,xmpio_single,fmarker,ierr)

   ! read nqlwl from the 2d record.
   positloc  = offset + bsize_frm + 9*xmpi_bsize_int
   call MPI_FILE_READ_AT(mpio_fh,positloc,nqlwl,1,MPI_INTEGER,statux,ierr)
   call wrtout(std_out, sjoin("nqlwl = ",itoa(nqlwl)))

   do isk=1,5
     call xmpio_read_frm(mpio_fh,offset,xmpio_single,fmarker,ierr)
   end do

   if (nqlwl>0) then  ! skip qlwl
     call xmpio_read_frm(mpio_fh,offset,xmpio_single,fmarker,ierr)
   end if

 case default
   MSG_BUG(sjoin('Wrong fform read:', itoa(fform)))
 end select

#else
 MSG_ERROR("hscr_mpio_skip cannot be used when MPI-IO is not enabled")
#endif

end subroutine hscr_mpio_skip
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/ioscr_qmerge
!! NAME
!! ioscr_qmerge
!!
!! FUNCTION
!!  Produce new file by merging the q-points stored in other files.
!!  This routine should be called by a single MPI process.
!!
!! INPUTS
!!  nfiles=Number of files to be merged.
!!  filenames(nfiles)=Paths of files to be merged.
!!  hscr_files(nfiles)<hscr_t>=Heades of the files to be merged.
!!  fname_out=Name of the file to be produced.
!!
!! OUTPUT
!!  ohscr<hscr_t>=The header of the output file.
!!
!! PARENTS
!!      mrgscr
!!
!! CHILDREN
!!      hscr_copy,hscr_io,hscr_print,read_screening,write_screening,wrtout
!!
!! SOURCE

subroutine ioscr_qmerge(nfiles, filenames, hscr_files, fname_out, ohscr)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfiles
 character(len=*),intent(in) :: fname_out
 type(hscr_t),intent(out) :: ohscr
!arrays
 character(len=*),intent(in) :: filenames(nfiles)
 type(hscr_t),intent(in) :: hscr_files(nfiles)

!Local variables-------------------------------
!scalars
 integer,parameter :: rdwr2=2,master=0
 integer :: iqibz,ifound,ifile,iqf,ount,iomode,fform_merge,comm,nomega4m,npwe4m,iqiA,ierr
 character(len=500) :: msg
 character(len=nctk_slen) :: varname
 type(abifile_t) :: abifile
!arrays
 integer,allocatable :: merge_table(:,:)
 real(dp) :: qdiff(3)
 complex(gwpc),allocatable :: epsm1(:,:,:,:)

! *************************************************************************

 comm = xmpi_comm_self

 if (file_exists(fname_out)) then
   MSG_ERROR(sjoin("Cannot overwrite existing file:", fname_out))
 end if

 ! Merge the headers creating the full list of q-points.
 call hscr_merge(Hscr_files(1:nfiles), ohscr)
 call hscr_print(ohscr, header='Header of the final file', unit=std_out, prtvol=1)

 ! For each q to be merged, save the index of the file where q is stored as well as its sequential index.
 ! Useful to do the merge point-by-point thus avoiding the allocation of the entire epsm1 array.
 ABI_MALLOC(merge_table,(ohscr%nqibz,2))
 do iqibz=1,ohscr%nqibz
   ifound=0
   fl: do ifile=1,nfiles
     do iqf=1,Hscr_files(ifile)%nqibz
       qdiff(:)=ohscr%qibz(:,iqibz)-Hscr_files(ifile)%qibz(:,iqf)
       if (all(abs(qdiff) < GW_TOLQ)) then
         merge_table(iqibz,1)=ifile
         merge_table(iqibz,2)=iqf
         ifound=ifound+1
         write(msg,'(a,3f12.6,2a)')'. q-point:',ohscr%qibz(:,iqibz),' will be taken from ',TRIM(filenames(ifile))
         call wrtout(std_out,msg,'COLL')
         EXIT fl
       end if
     end do
   end do fl
   ! Check if q-point has been found, multiple q-points not allowed.
   ABI_CHECK(ifound == 1, 'ifound/=1')
 end do

 iomode = IO_MODE_FORTRAN; if (endswith(fname_out, ".nc")) iomode = IO_MODE_ETSF
 if (iomode == IO_MODE_FORTRAN) then
   if (open_file(fname_out,msg,newunit=ount,status='new',form='unformatted') /= 0) then
     MSG_ERROR(msg)
   end if
 else
#ifdef HAVE_NETCDF
   NCF_CHECK(nctk_open_create(ount, fname_out, comm))
#endif
 end if

 ! Write the header.
 fform_merge = hscr_files(1)%fform
 abifile = abifile_from_fform(fform_merge)
 if (abifile%fform == 0) then
   MSG_ERROR(sjoin("Cannot find any abifile object associated to fform:", itoa(fform_merge)))
 end if
 varname = abifile%varname

 if (any(hscr_files(:)%fform /= hscr_files(1)%fform)) then
   write(std_out,*)"fforms: ",hscr_files(:)%fform
   MSG_ERROR("Files to be merged have different fform. Cannot merge data")
 end if

 call hscr_io(ohscr,fform_merge,rdwr2,ount,comm,master,iomode)

 npwe4m   = ohscr%npwe
 nomega4m = ohscr%nomega

 ABI_MALLOC_OR_DIE(epsm1,(npwe4m,npwe4m,nomega4m,1), ierr)

 do iqibz=1,ohscr%nqibz
   ifile = merge_table(iqibz,1)
   iqiA  = merge_table(iqibz,2)
   call read_screening(varname,filenames(ifile),npwe4m,1,nomega4m,epsm1,iomode,comm,iqiA=iqiA)
   call write_screening(varname,ount,iomode,npwe4m,nomega4m,iqibz,epsm1)
 end do

 ABI_FREE(epsm1)
 ABI_FREE(merge_table)

 if (iomode == IO_MODE_FORTRAN) then
   close(ount)
 else
#ifdef HAVE_NETCDF
   NCF_CHECK(nf90_close(ount))
#endif
 end if

 write(msg,'(3a)')ch10,' ==== Files have been merged successfully === ',ch10
 call wrtout(std_out,msg,'COLL')

end subroutine ioscr_qmerge
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/ioscr_qrecover
!! NAME
!! ioscr_qrecover
!!
!! FUNCTION
!!  Recover q-points from a corrupted file produced e.g. from an interrupted run
!!  This routine should be called by a single MPI process.
!!
!! INPUTS
!!  path=Corrupted file.
!!  nqrec=Number of q-points to recover.
!!  fname_out=Name of the file to be produced.
!!
!! OUTPUT
!!  Output is written to file.
!!
!! PARENTS
!!      mrgscr
!!
!! CHILDREN
!!      hscr_copy,hscr_io,hscr_print,read_screening,write_screening,wrtout
!!
!! SOURCE

subroutine ioscr_qrecover(ipath, nqrec, fname_out)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqrec
 character(len=*),intent(in) :: ipath,fname_out

!Local variables-------------------------------
!scalars
 integer,parameter :: rdwr2=2,master=0
 integer :: iqiA,nqibzA,nomega_asked,unt,npwe_asked,iomode,comm,fform1,ifform,ierr
 character(len=500) :: msg
 character(len=nctk_slen) :: varname
 type(hscr_t) :: hscr_recov,hscr
 type(abifile_t) :: abifile
!arrays
 complex(gwpc),allocatable :: epsm1(:,:,:,:)

! *************************************************************************

 comm = xmpi_comm_self

 call wrtout(std_out, sjoin(". Recovering q-points in file:", ipath))
 call wrtout(std_out, sjoin(". Data written to file:", fname_out))

 if (file_exists(fname_out)) then
   MSG_ERROR(sjoin("Cannot overwrite existing file:", fname_out))
 end if

 ! Find iomode from file extension and open output file.
 if (endswith(fname_out, ".nc")) then
   iomode = IO_MODE_ETSF
#ifdef HAVE_NETCDF
   NCF_CHECK(nctk_open_create(unt, fname_out, comm))
#else
   MSG_ERROR("Netcdf support is not available")
#endif
 else
   iomode = IO_MODE_FORTRAN
   if (open_file(fname_out, msg, newunit=unt, status='new', form='unformatted') /= 0) then
     MSG_ERROR(msg)
   end if
 end if

 ! Read header.
 call hscr_from_file(hscr, ipath, ifform, comm)
 ABI_CHECK(ifform /= 0, sjoin("fform = 0 while reading:", ipath))

 if (nqrec < 1 .or. nqrec > hscr%nqibz) then
   MSG_ERROR(sjoin("Wrong input. nqibz on file:", itoa(hscr%nqibz)))
 end if

 ! Copy header
 call hscr_copy(hscr, hscr_recov)

 ! Change dimensions and arrays associated to nqibz.
 hscr_recov%nqibz = nqrec
 ABI_FREE(hscr_recov%qibz)
 ABI_MALLOC(hscr_recov%qibz, (3,nqrec))
 hscr_recov%qibz = hscr%qibz(:,1:nqrec)

 call hscr_print(hscr_recov,header="Header of the new SCR file",unit=std_out,prtvol=1)

 ! Write the header of the recovered file.
 fform1 = hscr%fform

 abifile = abifile_from_fform(fform1)
 if (abifile%fform == 0) then
    MSG_ERROR(sjoin("Cannot find any abifile object associated to fform:", itoa(fform1)))
 end if
 varname = abifile%varname

 call hscr_io(hscr_recov,fform1,rdwr2,unt,comm,master,iomode)

 nqibzA=1; nomega_asked=hscr%nomega; npwe_asked=hscr%npwe

 ABI_MALLOC_OR_DIE(epsm1,(npwe_asked,npwe_asked,nomega_asked,1), ierr)

 do iqiA=1,hscr_recov%nqibz
   call read_screening(varname,ipath,npwe_asked,nqibzA,nomega_asked,epsm1,iomode,comm,iqiA=iqiA)
   call write_screening(varname,unt,iomode,npwe_asked,nomega_asked,iqiA,epsm1)
 end do

 if (iomode == IO_MODE_FORTRAN) close(unt)
#ifdef HAVE_NETCDF
 if (iomode == IO_MODE_ETSF) then
   NCF_CHECK(nf90_close(unt))
 end if
#endif

 ABI_FREE(epsm1)
 call hscr_free(hscr)
 call hscr_free(hscr_recov)

 call wrtout(std_out,"Recovery completed",'COLL')

end subroutine ioscr_qrecover
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/ioscr_wmerge
!! NAME
!! ioscr_wmerge
!!
!! FUNCTION
!!  Produce new file by merging the frequencies stored in other files.
!!  This routine should be called by a single MPI process.
!!
!! INPUTS
!!  nfiles=Number of files to be merged.
!!  filenames(nfiles)=Paths of files to be merged.
!!  hscr_files(nfiles)<hscr_t>=Heades of the files to be merged.
!!  fname_out=Name of the file to be produced.
!!
!! OUTPUT
!!  ohscr<hscr_t>=The header of the output file.
!!
!! PARENTS
!!      mrgscr
!!
!! CHILDREN
!!      hscr_copy,hscr_io,hscr_print,read_screening,write_screening,wrtout
!!
!! SOURCE

subroutine ioscr_wmerge(nfiles, filenames, hscr_file, freqremax, fname_out, ohscr)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfiles
 real(dp),intent(in) :: freqremax
 character(len=*),intent(in) :: fname_out
 type(hscr_t),intent(out) :: ohscr
!arrays
 character(len=*),intent(in) :: filenames(nfiles)
 type(hscr_t),intent(in) :: hscr_file(nfiles)

!Local variables-------------------------------
!scalars
 integer,parameter :: rdwr2=2,master=0
 integer :: ii,iqibz,ifile,ount,iomode,fform_merge,comm,nomega4m
 integer :: nfreq_tot,nfreqre,nfreqim,ifrq,npwe4mI,npwe4mJ,ierr
 character(len=500) :: msg
 logical :: skip
 character(len=nctk_slen) :: varname
 type(abifile_t) :: abifile
!arrays
 integer,allocatable :: freq_indx(:,:),ifile_indx(:),pos_indx(:),i_temp(:),i2_temp(:,:)
 real(dp),allocatable :: real_omega(:),imag_omega(:)
 complex(gwpc),allocatable :: epsm1(:,:,:,:),epsm1_temp(:,:,:,:)
 complex(dpc),allocatable :: omega_storage(:)

! *************************************************************************

 comm = xmpi_comm_self

 ! Check that q-point sets are the same
 do ifile=1,nfiles
   do ii=1,nfiles
     if (ii == ifile) CYCLE
     if (Hscr_file(ifile)%nqibz /= Hscr_file(ii)%nqibz) then
       MSG_ERROR(' One or more files do not have the same number of q-points!')
     end if
     do iqibz=1,Hscr_file(1)%nqibz
       if (ABS(SUM(Hscr_file(ifile)%qibz(:,iqibz)-Hscr_file(ii)%qibz(:,iqibz))) > tol6) then
         MSG_ERROR('Q-point set differs between one or more files!')
       end if
     end do
   end do
 end do

 ! nfreq_tot here is the total *possible* number of freq.
 nfreq_tot=0
 do ifile=1,nfiles
   nfreq_tot = nfreq_tot + Hscr_file(ifile)%nomega
 end do

 ABI_MALLOC(omega_storage, (nfreq_tot))
 ABI_MALLOC(freq_indx, (nfreq_tot,nfiles))
 ABI_MALLOC(ifile_indx, (nfreq_tot))
 omega_storage = CMPLX(-one,-one); freq_indx = 0; ifile_indx = 0

 ! Calculate the total number of real freq and store
 nfreqre = 0
 do ifile=1,nfiles
   do ifrq=1,Hscr_file(ifile)%nomega
     skip = .FALSE.
     ! Check whether to skip this point
     if (AIMAG(Hscr_file(ifile)%omega(ifrq)) > tol16) skip = .TRUE.
     if (REAL(Hscr_file(ifile)%omega(ifrq)) > freqremax) skip = .TRUE.
     ! Check for repetition or non-monotonic points
     if (nfreqre>1) then
       do ii=1,nfreqre
         if (ABS(REAL(Hscr_file(ifile)%omega(ifrq)) - REAL(omega_storage(ii))) < tol6) skip = .TRUE.
       end do
     end if
     if (skip) CYCLE
     nfreqre = nfreqre + 1
     ! Store (complex) frequency and index
     omega_storage(nfreqre) = Hscr_file(ifile)%omega(ifrq)
     ifile_indx(nfreqre) = ifile
     freq_indx(nfreqre,ifile) = ifrq
     write(std_out,'(a,es16.6,a,i0,2(a,i0))')&
       ' Found real frequency: ',REAL(omega_storage(nfreqre))*Ha_eV,' [eV], number: ',nfreqre,&
       ', in file: ',ifile,' local index: ',ifrq
   end do
 end do

 if (nfreqre > 0) then
   ! Sort real frequencies
   ABI_MALLOC(real_omega,(nfreqre))
   ABI_MALLOC(pos_indx,(nfreqre))
   ABI_MALLOC(i_temp,(nfreqre))
   ABI_MALLOC(i2_temp,(nfreqre,nfiles))
   real_omega(1:nfreqre) = REAL(omega_storage(1:nfreqre)) ! Copy real frequencies to temp. sorting array
   do ii=1,nfreqre ! Set up indexing array
     pos_indx(ii) = ii
   end do

   ! Sort frequencies while keeping track of index
   call sort_dp(nfreqre,real_omega,pos_indx,tol16)
   i_temp(1:nfreqre) = ifile_indx(1:nfreqre)
   i2_temp(1:nfreqre,1:nfiles) = freq_indx(1:nfreqre,1:nfiles)

   ! Copy sorted frequencies plus file and frequency index
   do ii=1,nfreqre
     omega_storage(ii) = CMPLX(real_omega(ii),zero)
     ifile_indx(ii) = i_temp(pos_indx(ii))
     freq_indx(ii,1:nfiles) = i2_temp(pos_indx(ii),1:nfiles)
   end do

   ABI_FREE(real_omega)
   ABI_FREE(pos_indx)
   ABI_FREE(i_temp)
   ABI_FREE(i2_temp)
 end if

 ! Check imaginary frequencies and store them
 nfreqim = 0
 do ifile=1,nfiles
   do ifrq=1,Hscr_file(ifile)%nomega
     if (REAL(Hscr_file(ifile)%omega(ifrq)) > tol8) CYCLE
     if (AIMAG(Hscr_file(ifile)%omega(ifrq)) < tol8) CYCLE
     nfreqim = nfreqim + 1
     omega_storage(nfreqre+nfreqim) = Hscr_file(ifile)%omega(ifrq)
     ifile_indx(nfreqre+nfreqim) = ifile
     freq_indx(nfreqre+nfreqim,ifile) = ifrq
     write(std_out,'(a,es16.6,a,i0,2(a,i0))')&
      ' Found imaginary frequency: ',AIMAG(omega_storage(nfreqre+nfreqim))*Ha_eV,' [eV], number: ',nfreqim,&
      ', in file: ',ifile,' local index: ',ifrq
   end do
 end do

 ! Sort imaginary frequencies
 ABI_MALLOC(imag_omega,(nfreqim))
 ABI_MALLOC(pos_indx,(nfreqim))
 ABI_MALLOC(i_temp,(nfreqim))
 ABI_MALLOC(i2_temp,(nfreqim,nfiles))

 ! Copy imaginary frequencies to temp. sorting array
 imag_omega(1:nfreqim) = AIMAG(omega_storage(nfreqre+1:nfreqre+nfreqim))
 do ii=1,nfreqim ! Set up indexing array
   pos_indx(ii) = ii
 end do

 ! Sort frequencies while keeping track of index
 call sort_dp(nfreqim,imag_omega,pos_indx,tol16)
 i_temp(1:nfreqim) = ifile_indx(nfreqre+1:nfreqre+nfreqim)
 i2_temp(1:nfreqim,1:nfiles) = freq_indx(nfreqre+1:nfreqre+nfreqim,1:nfiles)

 ! Copy sorted frequencies plus file and frequency index
 do ii=1,nfreqim
   omega_storage(nfreqre+ii) = CMPLX(zero,imag_omega(ii))
   ifile_indx(nfreqre+ii) = i_temp(pos_indx(ii))
   freq_indx(nfreqre+ii,1:nfiles) = i2_temp(pos_indx(ii),1:nfiles)
 end do

 ABI_FREE(imag_omega)
 ABI_FREE(pos_indx)
 ABI_FREE(i_temp)
 ABI_FREE(i2_temp)

 nfreq_tot = nfreqre + nfreqim ! Here nfreq_tot becomes the *true* number of freq
 write(std_out,'(2a,i0,a)') ch10,' Merging ',nfreq_tot,' frequencies.'
 write(std_out,'(2(a,i0),2a)') ' ',nfreqre,' real, and ',nfreqim,' imaginary.',ch10

 ! Copy old header
 call hscr_copy(Hscr_file(1),ohscr)

 ! TODO: hscr_wmerge
 ! Then modify entries for new frequency grid
 ohscr%nomega = nfreq_tot
 ABI_FREE(ohscr%omega)
 ABI_MALLOC(ohscr%omega,(nfreq_tot))
 ohscr%omega(1:nfreq_tot) = omega_storage(1:nfreq_tot)

 npwe4mI = ohscr%npwe*ohscr%nI
 npwe4mJ = ohscr%npwe*ohscr%nJ

 ! Print new header for info
 call hscr_print(ohscr,header='Header of the final file',unit=std_out,prtvol=1)

 if (file_exists(fname_out)) then
   MSG_ERROR(sjoin("Cannot overwrite existing file:", fname_out))
 end if

 if (endswith(fname_out, ".nc")) then
   iomode = IO_MODE_ETSF
#ifdef HAVE_NETCDF
   NCF_CHECK(nctk_open_create(ount, fname_out, comm))
#else
   MSG_ERROR("netcdf support is not activated")
#endif
 else
   iomode = IO_MODE_FORTRAN
   if (open_file(fname_out, msg, newunit=ount, status='new',form='unformatted') /= 0) then
     MSG_ERROR(msg)
   end if
 end if

 ! * Write the header.
 fform_merge = ohscr%fform

 abifile = abifile_from_fform(fform_merge)
 if (abifile%fform == 0) then
    MSG_ERROR(sjoin("Cannot find any abifile object associated to fform:", itoa(fform_merge)))
 end if
 varname = abifile%varname

 call hscr_io(ohscr,fform_merge,rdwr2,ount,comm,master,iomode)

 npwe4mI = ohscr%npwe*ohscr%nI
 npwe4mJ = ohscr%npwe*ohscr%nJ
 nomega4m = ohscr%nomega
 ABI_MALLOC_OR_DIE(epsm1,(npwe4mI,npwe4mJ,nomega4m,1), ierr)

 do iqibz=1,ohscr%nqibz

   do ifile=1,nfiles
     ! allocate temporary array
     npwe4mI = Hscr_file(ifile)%npwe*Hscr_file(ifile)%nI
     npwe4mJ = Hscr_file(ifile)%npwe*Hscr_file(ifile)%nJ
     nomega4m = Hscr_file(ifile)%nomega
     ABI_MALLOC_OR_DIE(epsm1_temp,(npwe4mI,npwe4mJ,nomega4m,1), ierr)

     ! read screening
     call read_screening(varname,filenames(ifile),npwe4mI,1,nomega4m,epsm1_temp,iomode,comm,iqiA=iqibz)

     ! Copy matrices for relevant frequencies
     do ifrq=1,nfreq_tot
       if (ifile_indx(ifrq)==ifile) then
         epsm1(:,:,ifrq,1)=epsm1_temp(:,:,freq_indx(ifrq,ifile),1)
       end if
     end do

     !Add imaginary part if needed
     !if (indx_imfreq_file==ifile) then
     !do ifrq=nfreqre+1,ohscr%nomega
     !epsm1(:,:,ifrq,1)=epsm1_temp(:,:,freq_indx(ifrq,ifile),1)
     !end do
     !end if

     ABI_FREE(epsm1_temp)
   end do !ifile

   ! Write data.
   npwe4mI = ohscr%npwe*ohscr%nI
   nomega4m = ohscr%nomega
   call write_screening(varname,ount,iomode,npwe4mI,nomega4m,iqibz,epsm1)
 end do ! iqibz

 ABI_FREE(epsm1)
 ABI_FREE(omega_storage)
 ABI_FREE(freq_indx)
 ABI_FREE(ifile_indx)

 if (iomode == IO_MODE_FORTRAN) then
   close(ount)
 else
#ifdef HAVE_NETCDF
   NCF_CHECK(nf90_close(ount))
#endif
 end if

 write(msg,'(3a)')ch10,' ==== Files have been merged successfully === ',ch10
 call wrtout(std_out,msg,'COLL')

end subroutine ioscr_wmerge
!!***

!----------------------------------------------------------------------

!!****f* m_io_screening/ioscr_wremove
!! NAME
!! ioscr_wremove
!!
!! FUNCTION
!!  Produce new file by removing selected frequencies in the initial file `inpath`.
!!  This routine should be called by a single MPI process.
!!
!! INPUTS
!!  inpath=Input file
!!  ihscr<hscr_t>=Headerf of the input file.
!!  fname_out=Output file.
!!  nfreq_tot=Number of frequencies in new file.
!!  freq_indx(nfreq_tot)=Index of frequency to be kept in input file.
!!
!! OUTPUT
!!  ohscr<hscr_t>=The header of the output file.
!!
!! PARENTS
!!      mrgscr
!!
!! CHILDREN
!!      hscr_copy,hscr_io,hscr_print,read_screening,write_screening,wrtout
!!
!! SOURCE

subroutine ioscr_wremove(inpath, ihscr, fname_out, nfreq_tot, freq_indx, ohscr)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfreq_tot
 character(len=*),intent(in) :: inpath,fname_out
 type(hscr_t),intent(in) :: ihscr
 type(hscr_t),intent(out) :: ohscr
!arrays
 integer,intent(in) :: freq_indx(nfreq_tot)

!Local variables-------------------------------
!scalars
 integer,parameter :: rdwr2=2,master=0
 integer :: iqibz,fform_merge,comm,nomega4m,ierr
 integer :: ifrq,npwe4mI,npwe4mJ,iomode,ount
 character(len=500) :: msg
 character(len=nctk_slen) :: varname
 type(abifile_t) :: abifile
!arrays
 complex(gwpc),allocatable :: epsm1(:,:,:),epsm1_temp(:,:,:)

! *************************************************************************

 comm = xmpi_comm_self

 ! check ifreq_idx
 ABI_CHECK(nfreq_tot > 0, "nfreq_tot <= 0!")
 if (all(freq_indx == 0)) MSG_ERROR("all(freq_indx == 0)")

 ! Copy the old header
 call hscr_copy(ihscr, ohscr)

 ! Then modify entries for new frequency grid.
 ohscr%nomega = nfreq_tot
 ABI_FREE(ohscr%omega)
 ABI_MALLOC(ohscr%omega,(nfreq_tot))
 do ifrq=1,nfreq_tot
   ohscr%omega(ifrq) = ihscr%omega(freq_indx(ifrq))
 end do

 npwe4mI = ohscr%npwe*ohscr%nI
 npwe4mJ = ohscr%npwe*ohscr%nJ

 ! Print new header for info
 call hscr_print(ohscr,header='Header of the final file',unit=std_out,prtvol=1)

 ! Open output file.
 if (endswith(fname_out, ".nc")) then
   iomode = IO_MODE_ETSF
#ifdef HAVE_NETCDF
   NCF_CHECK(nctk_open_create(ount, fname_out, comm))
#else
   MSG_ERROR("Netcdf support is not available")
#endif
 else
   iomode = IO_MODE_FORTRAN
   if (open_file(fname_out, msg, newunit=ount, status='new', form='unformatted') /= 0) then
     MSG_ERROR(msg)
   end if
 end if

 ! Write the header.
 fform_merge = ohscr%fform
 abifile = abifile_from_fform(fform_merge)
 if (abifile%fform == 0) then
    MSG_ERROR(sjoin("Cannot find any abifile object associated to fform:", itoa(fform_merge)))
 end if
 varname = abifile%varname

 call hscr_io(ohscr,fform_merge,rdwr2,ount,comm,master,iomode)

 npwe4mI = ohscr%npwe*ohscr%nI; npwe4mJ = ohscr%npwe*ohscr%nJ
 nomega4m = ohscr%nomega

 ABI_MALLOC_OR_DIE(epsm1, (npwe4mI,npwe4mJ,nomega4m), ierr)

 do iqibz=1,ohscr%nqibz
   ! allocate temporary array
   npwe4mI = ihscr%npwe * ihscr%nI
   npwe4mJ = ihscr%npwe * ihscr%nJ
   nomega4m = ihscr%nomega
   ABI_MALLOC_OR_DIE(epsm1_temp,(npwe4mI,npwe4mJ,nomega4m), ierr)

   ! read full screening matrix for this q-point
   call read_screening(varname,inpath,npwe4mI,1,nomega4m,epsm1_temp,iomode,comm,iqiA=iqibz)

   ! Copy relevant frequencies
   do ifrq=1,nfreq_tot
     epsm1(:,:,ifrq) = epsm1_temp(:,:,freq_indx(ifrq))
   end do

   ABI_FREE(epsm1_temp)

   npwe4mI = ohscr%npwe*ohscr%nI; nomega4m = ohscr%nomega
   call write_screening(varname,ount,iomode,npwe4mI,nomega4m,iqibz,epsm1)
 end do ! iqibz

 ABI_FREE(epsm1)

 if (iomode == IO_MODE_FORTRAN) then
   close(ount)
 else
#ifdef HAVE_NETCDF
   NCF_CHECK(nf90_close(ount))
#endif
 end if

 write(msg,'(3a)')ch10,' ==== Frequencies have been removed successfully === ',ch10
 call wrtout(std_out,msg,'COLL')

end subroutine ioscr_wremove
!!***

END MODULE m_io_screening
!!***
