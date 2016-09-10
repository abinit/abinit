!{\src2tex{textfont=tt}}
!!****f* ABINIT/inpgkk
!! NAME
!! inpgkk
!!
!! FUNCTION
!! read in gkk file and return eigenvalue matrix
!! Only works for a single gkk matrix (1 perturbation and qpoint) in the file
!! like the files produced by outgkk
!!
!! COPYRIGHT
!! Copyright (C) 2009-2016 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  
!!  filegkk= filename
!!
!! OUTPUT
!!  eigen1 = response function 1st order eigenvalue matrix
!!
!! PARENTS
!!      read_el_veloc
!!
!! CHILDREN
!!      hdr_fort_read,hdr_free,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine inpgkk(eigen1,filegkk,hdr1)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_hdr

 use m_io_tools,        only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'inpgkk'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: filegkk
 type(hdr_type), intent(out) :: hdr1
!arrays
 real(dp),allocatable,intent(out) :: eigen1(:)

!Local variables-------------------------------
!scalars
 integer :: bantot1
 integer :: isppol, ikpt, mband, ikb
 integer :: unitgkk, fform, ierr, n1wf, i1wf
 type(hdr_type) :: hdr0
 real(dp), allocatable :: eigen(:)
 character(len=500) :: message

! *************************************************************************

 if (open_file(filegkk,message,newunit=unitgkk,form='unformatted',status='old') /= 0) then
   MSG_ERROR(message)
 end if


!read in header of GS file and eigenvalues
 call hdr_fort_read(hdr0, unitgkk, fform)
 ABI_CHECK(fform /= 0, "hdr_fort_read returned fform == 0")
 
 mband = maxval(hdr0%nband(:))
 ABI_ALLOCATE(eigen,(mband))
 call wrtout(std_out,'inpgkk : try to reread GS eigenvalues','COLL')

 do isppol=1,hdr0%nsppol
   do ikpt=1,hdr0%nkpt
     read (unitgkk,IOSTAT=ierr) eigen(1:hdr0%nband(ikpt))
     ABI_CHECK(ierr==0,'reading eigen from gkk file')
   end do
 end do

 read(unitgkk,IOSTAT=ierr) n1wf
 ABI_CHECK(ierr==0,"reading n1wf from gkk file")

 ABI_DEALLOCATE(eigen)
 call hdr_free(hdr0)

 if (n1wf > 1) then
   write(message,'(3a)')&
&   'several 1wf records were found in the file,',ch10, &
&   'which is not allowed for reading with this routine'
   MSG_ERROR(message)
 end if

!read in header of 1WF file
 call hdr_fort_read(hdr1, unitgkk, fform)
 if (fform == 0) then
   write(message,'(a,i0,a)')' 1WF header number ',i1wf,' was mis-read. fform == 0'
   MSG_ERROR(message)
 end if

 bantot1 = 2*hdr1%nsppol*hdr1%nkpt*mband**2
 ABI_ALLOCATE(eigen1, (bantot1))


!retrieve 1WF <psi_k+q | H | psi_k> from gkk file and echo to output
 ikb = 0
 do isppol=1,hdr1%nsppol
   do ikpt=1,hdr1%nkpt
     read (unitgkk,IOSTAT=ierr) eigen1(ikb+1:ikb+2*hdr1%nband(ikpt)**2)
     ikb = ikb + 2*hdr1%nband(ikpt)**2
     if (ierr /= 0) then
       write(message,'(a,2i0)')'reading eigen1 from gkk file, spin, kpt_idx',isppol,ikpt
       MSG_ERROR(message)
     end if
   end do
 end do

 close(unitgkk)

end subroutine inpgkk
!!***
