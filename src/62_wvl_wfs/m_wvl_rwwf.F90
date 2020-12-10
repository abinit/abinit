!!****m* ABINIT/m_wvl_rwwf
!! NAME
!!  m_wvl_rwwf
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DC)
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

module m_wvl_rwwf

 use defs_basis

 use defs_wvltypes
 use m_wffile
 use m_errors
 use m_abicore
 use m_hdr
 use m_xmpi
 use m_dtset

 use defs_abitypes,  only : MPI_type
 use m_geometry,     only : xred2xcart

 implicit none

 private
!!***

 public :: wvl_read
 public :: wvl_write
!!***

contains
!!***

!!****f* ABINIT/wvl_read
!! NAME
!! wvl_read
!!
!! FUNCTION
!! Simple wrapper around the read disk methods of BigDFT for wavefunctions.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=input variables.
!!  mpi_enreg=information about MPI parallelization
!!  option= -2 for write with BigDFT format,
!!          -1 for reading wavelets coefficients with BigDFT format,
!!          2 for write,
!!          1 for read.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  wff <type(wffile_type)>=struct info for wavefunction
!!  wfs <type(wvl_wf_type)>=wavefunctions information for wavelets.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!
!! PARENTS
!!      m_wvl_wfsinp
!!
!! CHILDREN
!!      etsf_io_basisdata_put,etsf_io_electrons_put,etsf_io_main_put
!!      writemywaves,writeonewave,wrtout,xred2xcart
!!
!! SOURCE

subroutine wvl_read(dtset, hdr0, hdr, mpi_enreg, option, rprimd, wff, wfs, wvl, xred)

#if defined HAVE_BIGDFT
  use BigDFT_API, only: readonewave, reformatonewave, readmywaves, &
&                       WF_FORMAT_NONE
#endif
  implicit none

!Arguments -------------------------------
  !scalars
  integer, intent(in)                       :: option
  type(dataset_type), intent(in)            :: dtset
  type(hdr_type), intent(in)                :: hdr0
  type(hdr_type), intent(in)                :: hdr
  type(MPI_type), intent(in)                :: mpi_enreg
  type(wffile_type),intent(in)              :: wff
  type(wvl_wf_type), intent(inout)          :: wfs
  type(wvl_internal_type), intent(in)       :: wvl
  !arrays
  real(dp), intent(in)                      :: rprimd(3, 3)
  real(dp), intent(in)                      :: xred(3, dtset%natom)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
  character(len = 500)  :: message
  integer               :: iBand, bandSize
  integer               :: comm,me
  real(dp), allocatable :: xcart(:,:), psifscf(:,:,:)
  real(dp), allocatable :: xcart_old(:,:)
#endif

! *********************************************************************

#if defined HAVE_BIGDFT

 if (abs(option) /= 1) then
   write(message,'(a,a,a,i0,a)')&
&   '  Option argument is wrong,', ch10, &
&   '  awaited values are -1 or  1 but option = ', option, '.'
   ABI_BUG(message)
 end if

 comm=mpi_enreg%comm_wvl
 me=xmpi_comm_rank(comm)
!Store xcart for each atom
 ABI_ALLOCATE(xcart,(3, dtset%natom))
 ABI_ALLOCATE(xcart_old,(3, dtset%natom))
 call xred2xcart(dtset%natom, rprimd, xcart, xred)

 write(message,'(2a)') ch10,' wvl_read:  read wavefunctions from file.'
 call wrtout(std_out,message,'COLL')

 if (option > 0) then
   bandSize = wfs%ks%lzd%Glr%wfd%nvctr_c + 7 * wfs%ks%lzd%Glr%wfd%nvctr_f
!  Read in the ABINIT way.
   if (wff%iomode == IO_MODE_FORTRAN .or. (wff%iomode == IO_MODE_FORTRAN_MASTER .and. wff%master==wff%me)) then
     ABI_ALLOCATE(psifscf,(wvl%Glr%d%n1i, wvl%Glr%d%n2i, wvl%Glr%d%n3i))
     do iBand = 1, dtset%mband * dtset%nsppol, 1
       call readonewave(wff%unwff, .false., iBand, me, &
&       wvl%Glr%d%n1, wvl%Glr%d%n2, wvl%Glr%d%n3, &
&       wvl%h(1), wvl%h(2), wvl%h(3), wvl%atoms, &
&       wfs%ks%lzd%Glr%wfd, xcart_old, xcart, &
&       wfs%ks%psi(bandSize * (iBand - me * wfs%ks%orbs%norbp - 1) + 1: &
&       bandSize * (iBand - me * wfs%ks%orbs%norbp - 1) + bandSize), &
&       wfs%ks%orbs%eval(iBand), psifscf)
     end do
     ABI_DEALLOCATE(psifscf)

   else
     write(message,'(4a,i0,a)') ch10,&
&     '  wff%iomode argument is wrong,', ch10, &
&     '  awaited values are -1, 0 (or 3 if netcdf/etsf_io is available) but value = ', wff%iomode, '.'
     ABI_BUG(message)
   end if
 else
   call readmywaves(me, "wavefunctions", WF_FORMAT_NONE, wfs%ks%orbs, &
&   wvl%Glr%d%n1, wvl%Glr%d%n2, wvl%Glr%d%n3, &
&   wvl%h(1), wvl%h(2), wvl%h(3), wvl%atoms, &
&   xcart_old, xcart, wfs%ks%lzd%Glr%wfd, wfs%ks%psi)
 end if

 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(xcart_old)
#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) option,dtset%nstep,hdr0%ecut,hdr%ecut,mpi_enreg%nproc,wff%me,&
& wfs%ks,wvl%h(1),rprimd(1,1),xred(1,1)
#endif

end subroutine wvl_read
!!***

!!****f* ABINIT/wvl_write
!! NAME
!! wvl_write
!!
!! FUNCTION
!! Simple wrapper around the write disk methods of BigDFT for wavefunctions.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=input variables.
!!  mpi_enreg=information about MPI parallelization
!!  option= -2 for write with BigDFT format,
!!          -1 for reading wavelets coefficients with BigDFT format,
!!          2 for write,
!!          1 for read.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  wff <type(wffile_type)>=struct info for wavefunction
!!  wfs <type(wvl_wf_type)>=wavefunctions information for wavelets.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! PARENTS
!!      m_iowf
!!
!! CHILDREN
!!      etsf_io_basisdata_put,etsf_io_electrons_put,etsf_io_main_put
!!      writemywaves,writeonewave,wrtout,xred2xcart
!!
!! SOURCE

subroutine wvl_write(dtset, eigen, mpi_enreg, option, rprimd, wff, wfs, wvl, xred)

#if defined HAVE_BIGDFT
  use BigDFT_API, only : writeonewave,writemywaves,WF_FORMAT_NONE
#endif
  implicit none

!Arguments -------------------------------
  !scalars
  integer, intent(in)                       :: option
  type(dataset_type), intent(in)            :: dtset
  type(MPI_type), intent(in)                :: mpi_enreg
  type(wffile_type),intent(in)              :: wff
  type(wvl_wf_type), intent(in)             :: wfs
  type(wvl_internal_type), intent(in)       :: wvl
  !arrays
  real(dp), intent(in), target              :: eigen(dtset%mband)
  real(dp), intent(in)                      :: rprimd(3, 3)
  real(dp), intent(in)                      :: xred(3, dtset%natom)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
  character(len = 500)  :: message
  integer               :: comm,me
  integer               :: iorb
  integer               :: iseg, nseg, ipsi, npsi
  real(dp), allocatable :: xcart(:,:)
#endif

! *********************************************************************

#if defined HAVE_BIGDFT

 if (abs(option) /= 2) then
   write(message,'(a,a,a,i0,a)')&
&   '  Option argument is wrong,', ch10, &
&   '  awaited values are -2 or  2 but option = ', option, '.'
   ABI_BUG(message)
 end if

 comm=mpi_enreg%comm_wvl
 me=xmpi_comm_rank(comm)
!Store xcart for each atom
 ABI_ALLOCATE(xcart,(3, dtset%natom))
 call xred2xcart(dtset%natom, rprimd, xcart, xred)

 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_write:  Write wavefunctions to file.'
 call wrtout(std_out,message,'COLL')

 if (option > 0) then
!  Write in the ABINIT way.
   if (wff%iomode == IO_MODE_FORTRAN .or. (wff%iomode == IO_MODE_FORTRAN_MASTER .and. wff%master==wff%me)) then
     iseg = wfs%ks%lzd%Glr%wfd%nseg_c
     nseg = wfs%ks%lzd%Glr%wfd%nseg_c + wfs%ks%lzd%Glr%wfd%nseg_f
     ipsi = wfs%ks%lzd%Glr%wfd%nvctr_c
     npsi = wfs%ks%lzd%Glr%wfd%nvctr_c + 7 * wfs%ks%lzd%Glr%wfd%nvctr_f
     do iorb = 1, dtset%mband
       call writeonewave(wff%unwff, .false., iorb, wvl%Glr%d%n1, &
&       wvl%Glr%d%n2, wvl%Glr%d%n3, &
&       wvl%h(1), wvl%h(2), wvl%h(3), dtset%natom, &
&       xcart, wfs%ks%lzd%Glr%wfd%nseg_c, wfs%ks%lzd%Glr%wfd%nvctr_c, &
&       wfs%ks%lzd%Glr%wfd%keygloc(:,1:iseg), wfs%ks%lzd%Glr%wfd%keyvloc(1:iseg), wfs%ks%lzd%Glr%wfd%nseg_f, &
&       wfs%ks%lzd%Glr%wfd%nvctr_f, wfs%ks%lzd%Glr%wfd%keygloc(:, iseg + 1:nseg), &
&       wfs%ks%lzd%Glr%wfd%keyvloc(iseg + 1:nseg), &
&       wfs%ks%psi(npsi * (iorb - me * wfs%ks%orbs%norbp - 1) + 1: &
&       npsi * (iorb - me * wfs%ks%orbs%norbp - 1) + ipsi), &
&       wfs%ks%psi(npsi * (iorb - me * wfs%ks%orbs%norbp - 1) + ipsi + 1: &
&       npsi * (iorb - me * wfs%ks%orbs%norbp - 1) + npsi), &
&       wfs%ks%orbs%eval(iorb))
     end do
   else
     write(message,'(3a,i0,a)')&
&     '  wff%iomode argument is wrong,', ch10, &
&     '  awaited values are -1, 0 (or 3 if netcdf/etsf_io is available) but value = ', wff%iomode, '.'
     ABI_BUG(message)
   end if
 else
!  Write in the BigDFT way.
   call  writemywaves(me, "wavefunctions", WF_FORMAT_NONE, wfs%ks%orbs, &
&   wvl%Glr%d%n1, wvl%Glr%d%n2, wvl%Glr%d%n3, &
&   wvl%h(1), wvl%h(2), wvl%h(3),wvl%atoms, &
&   xcart, wfs%ks%lzd%Glr%wfd, wfs%ks%psi)
 end if

 ABI_DEALLOCATE(xcart)

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) option,dtset%nstep,mpi_enreg%nproc,wff%me,&
& wfs%ks,wvl%h(1),eigen(1),rprimd(1,1),xred(1,1)
#endif

end subroutine wvl_write
!!***

end module m_wvl_rwwf
!!***
