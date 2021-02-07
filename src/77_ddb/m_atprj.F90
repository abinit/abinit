!!****m* ABINIT/m_atprj
!!
!! NAME
!! m_atprj
!!
!! FUNCTION
!! Module to output atomic projections of phonon modes
!!
!! COPYRIGHT
!! Copyright (C) 2011-2020 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_atprj

 use defs_basis
 use m_abicore
 use m_errors

 use m_io_tools, only : get_unit, open_file
 use m_fstrings, only : int2char4

 implicit none

 private
!!***

!!****t* m_atprj/atprj_type
!! NAME
!! atprj_type
!!
!! FUNCTION
!! Container for atomic projection file data
!!
!! SOURCE

type, public :: atprj_type

  integer :: natprj_bs
  integer :: natom

  integer, allocatable :: iatprj_bs(:)
  character(len=fnlen), allocatable :: filename(:,:)

end type atprj_type

public :: atprj_init
public :: atprj_print
public :: atprj_destroy

contains
!!***

!!****f* m_atprj/atprj_init
!!
!! NAME
!! atprj_init
!!
!! FUNCTION
!! initialize atprj datastructure
!!
!! INPUT
!! natom = number of atoms
!! natprj_bs = number of atoms to project on
!! iatprj_bs = indices of atoms to project on
!! outfile_radix = base file name for output files
!!
!! OUTPUT
!! t_atprj = container object for atomic projections
!!
!! PARENTS
!!      m_phonons
!!
!! CHILDREN
!!
!! SOURCE

subroutine atprj_init(t_atprj, natom, natprj_bs, iatprj_bs, outfile_radix)

 type(atprj_type), intent(out) :: t_atprj
 integer, intent(in) :: natom
 integer, intent(in) :: natprj_bs
 integer, intent(in) :: iatprj_bs(natprj_bs)
 character(len=*), intent(in) :: outfile_radix

!Local variables-------------------------------
!scalars
 integer :: iatom, imode, iunit
 character(len=10) :: imodestring, iatomstring
 character(len=500) :: msg
! *************************************************************************

 t_atprj%natprj_bs = natprj_bs
 t_atprj%natom = natom

 ABI_MALLOC(t_atprj%iatprj_bs,(natprj_bs))
 t_atprj%iatprj_bs = iatprj_bs

! for each phonon mode and atom for projection, open a file
 ABI_MALLOC(t_atprj%filename ,(3*natom,natprj_bs))
 iunit = get_unit()
 do imode = 1, 3*natom
   call int2char4(imode, imodestring)
   ABI_CHECK((imodestring(1:1)/='#'),'Bug: string length too short!')
   do iatom = 1, natprj_bs
     call int2char4(iatom, iatomstring)
     ABI_CHECK((iatomstring(1:1)/='#'),'Bug: string length too short!')
     t_atprj%filename(imode,iatom) = trim(outfile_radix)//"_mod"//trim(imodestring)//"_iat"//trim(iatomstring)
     if (open_file(t_atprj%filename(imode,iatom), msg, newunit=iunit, form="formatted", action="write") /= 0) then
       ABI_ERROR(msg)
     end if
     ! print header
     write (unit=iunit, fmt='(a)') '##'
     write (unit=iunit, fmt='(a,I6,a)') '##  This file contains abinit phonon frequencies for mode number ', &
&          imode, ' along a path in reciprocal space,'
     write (unit=iunit, fmt='(a,I6)') '##  the third column is the projection along atom number ',&
&          t_atprj%iatprj_bs(iatom)
     write (unit=iunit, fmt='(a)') '##'

     close (iunit)
   end do
 end do

end subroutine atprj_init
!!***

!!****f* m_atprj/atprj_print
!!
!! NAME
!! atprj_print
!!
!! FUNCTION
!! print out 1 line per atomic projection file
!!
!! INPUT
!! t_atprj = container object for atomic projections
!! phfrq = phonon frequencies for present q point
!! eigvec = eigenvectors for present q point
!!
!! OUTPUT
!!  writes to files
!!
!! PARENTS
!!      m_phonons
!!
!! CHILDREN
!!
!! SOURCE

subroutine atprj_print(t_atprj, iq, phfrq, eigvec)

!arguments
 integer, intent(in) :: iq
 type(atprj_type), intent(in) :: t_atprj
 real(dp), intent(in) :: phfrq(3*t_atprj%natom)
 real(dp), intent(in) :: eigvec(2,3,t_atprj%natom,3,t_atprj%natom)

!local variables
 integer :: jatom, idir, iatom, imode, iunit, jdir
 real(dp) :: proj

! *************************************************************************

 iunit = get_unit()
 do iatom = 1, t_atprj%natom
   do idir = 1, 3
     imode = idir + 3*(iatom-1)
     do jatom = 1, t_atprj%natprj_bs
       open (unit=iunit, file=t_atprj%filename(imode,jatom), position='append')
       write (unit=iunit, fmt='(a,I4,a)') '# atom ', jatom, ' sum of directions'
       proj = sum(eigvec(:,:,jatom,idir,iatom)**2)
       write (unit=iunit, fmt='(I6,2E20.10)') iq, phfrq(imode), proj

       do jdir = 1, 3
         write (unit=iunit, fmt='(2a,I4,a,I4)') ch10, '# atom ', jatom, ' directions ', jdir
         proj = sum(eigvec(:,jdir,jatom,idir,iatom)**2)
         write (unit=iunit, fmt='(I6,2E20.10)') iq, phfrq(imode), proj
       end do
       close (iunit)
     end do
   end do
 end do

end subroutine atprj_print
!!***

!!****f* m_atprj/atprj_destroy
!!
!! NAME
!! atprj_destroy
!!
!! FUNCTION
!! deallocate atomic projection datastructure and close files
!!
!! INPUT
!!
!! OUTPUT
!! t_atprj = container object for atomic projections
!!
!! PARENTS
!!      m_phonons
!!
!! CHILDREN
!!
!! SOURCE

subroutine atprj_destroy(t_atprj)

 type(atprj_type), intent(inout) :: t_atprj

 if (allocated(t_atprj%iatprj_bs)) then
   ABI_FREE(t_atprj%iatprj_bs)
 end if

 if (allocated(t_atprj%filename)) then
   ABI_FREE(t_atprj%filename)
 end if

end subroutine atprj_destroy

end module m_atprj
!!***
