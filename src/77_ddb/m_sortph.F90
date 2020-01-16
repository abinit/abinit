!!****m* ABINIT/m_sorth_ph
!! NAME
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (MVer, FDortu, MVeithen)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_sortph

 use defs_basis
 use m_abicore
 use m_errors

 use m_io_tools,   only : open_file

 implicit none

 private

 complex(dpc),save,allocatable :: eigvecLast(:,:)

 public :: end_sortph
 public :: sortph

 ! Logical units used to write data.
 integer,private,save :: udispl=-1,ufreq=-1
!!***

contains

!!****f* m_sortph/end_sortph
!! NAME
!! end_sortph
!!
!! FUNCTION
!! Deallocate array for sortph
!!
!! INPUTS
!!
!! OUTPUT
!!  Only deallocation
!!
!! NOTES
!!
!! PARENTS
!!      harmonic_thermo,m_phonons
!!
!! CHILDREN
!!
!! SOURCE
subroutine end_sortph()

 if (allocated(eigvecLast))  then
   ABI_DEALLOCATE(eigvecLast)
 end if

 if (ufreq /= -1) then
   close(ufreq); ufreq = -1
 end if
 if (udispl /= -1) then
   close(udispl); udispl = -1
 end if

end subroutine end_sortph
!!***

!!****f* m_sortph/sortph
!! NAME
!! sortph
!!
!! FUNCTION
!! Sort the energies in order to have fine phonon dispersion curves
!! It is best not to include the gamma point in the list
!!
!! MODIFIED
!! Takeshi Nishimatsu
!!
!! INPUTS
!!  eigvec(2*3*natom*3*natom)= contain
!!  the eigenvectors of the dynamical matrix.
!!  displ(2*3*natom*3*natom)= contain
!!   the displacements of atoms in cartesian coordinates.
!!   The first index means either the real or the imaginary part,
!!   The second index runs on the direction and the atoms displaced
!!   The third index runs on the modes.
!!  filnam=name of output files
!!   hacmm1,hartev,harthz,xkb= different conversion factors
!!  natom= number of atom
!!  phfrq(3*natom)= phonon frequencies in Hartree
!!
!! OUTPUT
!!  (only writing ?)
!!
!! NOTES
!! Called by one processor only
!!
!! PARENTS
!!      harmonic_thermo,m_phonons
!!
!! CHILDREN
!!
!! SOURCE

subroutine sortph(eigvec,displ,filnam, natom,phfrq)

!Arguments -----------------------------------
!scalars
integer,intent(in) :: natom
character(len=*),intent(in) :: filnam
!arrays
real(dp),intent(in) :: eigvec(2,3,natom,3,natom)
real(dp),intent(in) :: displ(2*3*natom*3*natom)
real(dp),intent(in) :: phfrq(3*natom)

!Local variables-------------------------------
!scalars
integer :: iatom,imode,j,idir1,idir2,ipert1,ipert2,i1,i2
character(len=fnlen) :: file_displ,file_freq
character(len=20) :: fmt_phfrq
character(len=500) :: msg
!arrays
logical     ::               mask(3*natom)
real(dp)    ::           phfrqNew(3*natom)
complex(dpc) ::           displIn(3*natom,3*natom)
complex(dpc) ::           displNew(3*natom,3*natom)
complex(dpc) ::           eigvecIn(3*natom,3*natom)
complex(dpc) ::          eigvecNew(3*natom,3*natom)
complex(dpc) ::   transpose_eigvec(3*natom,3*natom)
real(dp)    ::     abs_similarity(3*natom,3*natom)  !|<displNew|displLast>|
! *********************************************************************

do ipert2=1,natom
  do idir2=1,3
    i2=idir2+(ipert2-1)*3
    do ipert1=1,natom
      do idir1=1,3
        i1=idir1+(ipert1-1)*3
        eigvecIn(i1,i2)=cmplx(eigvec(1,idir1,ipert1,idir2,ipert2),eigvec(2,idir1,ipert1,idir2,ipert2))
        displIn(i1,i2)=cmplx(displ(1+2*(i1-1)+2*3*natom*(i2-1)),displ(2+2*(i1-1)+2*3*natom*(i2-1)))
      end do
    end do
  end do
end do

 if(.not.allocated(eigvecLast)) then
   file_freq  = trim(filnam)//".freq" !---------------------------------------------------
   write(std_out,'(a,a)' )' sortph : opening file ',trim(file_freq)
   if (open_file(file_freq,msg,newunit=ufreq,STATUS='replace',ACTION='write') /= 0) then
     MSG_ERROR(msg)
   end if
   file_displ = trim(filnam)//".displ" !--------------------------------------------------
   write(std_out,'(a,a)' )' sortph : opening file ',trim(file_displ)
   if (open_file(file_displ,msg,newunit=udispl,STATUS='replace',ACTION='write') /= 0) then
     MSG_ERROR(msg)
   end if
   ABI_ALLOCATE(eigvecLast,(3*natom,3*natom))
   phfrqNew(:)   =  phfrq(:)
   displNew(:,:) =  displIn(:,:)
   eigvecNew(:,:) = eigvecIn(:,:)
 else
!  Avoid gfortran 4.2.1 bug, with which you CANNOT conjg(transpose(displ))
   transpose_eigvec = transpose(eigvecIn)
   abs_similarity = abs(matmul(conjg(transpose_eigvec),eigvecLast))
   mask(:) = .true.
   phfrqNew(:)   =  phfrq(:)
   displNew(:,:) =  displIn(:,:)
   eigvecNew(:,:) = eigvecIn(:,:)
 end if


!Write frequencies in a file
 write(fmt_phfrq,'(a,i3,a)') '(', 3*natom, 'e18.10)'
 write(ufreq,fmt_phfrq) (phfrqNew(j),j=1,3*natom)

!write displacements in a file
! NB: sqrt still returns a complex number could be using modulus or something simpler
   do imode=1,3*natom
     do iatom=1,natom
       write(udispl,'(e18.10)') &
       real(sqrt(displNew(3*(iatom-1)+1,imode) * conjg(displNew(3*(iatom-1)+1,imode)) + &
&      displNew(3*(iatom-1)+2,imode) * conjg(displNew(3*(iatom-1)+2,imode)) + &
&      displNew(3*(iatom-1)+3,imode) *  conjg(displNew(3*(iatom-1)+3,imode)) ))
     end do
   end do

 eigvecLast(:,:) = eigvecNew(:,:)

end subroutine sortph
!!***

end module m_sortph
!!***
