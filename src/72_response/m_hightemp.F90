!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_hightemp
!! NAME
!!  m_hightemp
!!
!! FUNCTION
!!  This module provides routines to run computations at very high temperatures
!!  with orbitals.
!!
!! COPYRIGHT
!!  Copyright (C) 2018-2019 ABINIT group (??)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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

module hightemp
  
  use defs_basis
  use m_io_tools
  use m_errors

  implicit none

  private
  
  public :: prt_eigocc

contains

!!****f* ABINIT/prt_eigocc
!! NAME
!! prt_eigocc
!!
!! FUNCTION
!! Printing in a _EIGOCC file many informations like Fermi energy, Average Vxc, Unit cell volume,
!! electronic temperature (tsmear when occopt=3), eigenvalues and occupation for each k-point,
!! wheight of the k-point, number of bands, number of k-points...
!! This file is intended to be used for custom DOS computation with external tools for example.
!!
!! COPYRIGHT
!! Copyright (C) 2018-2019 ABINIT group (??)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! TODO
!!
!! INPUTS
!! chksymbreak= if 1, will check whether the k point grid is symmetric, and stop if not.
!! gmet(3,3)=reciprocal space metric (bohr**-2).
!! iout=if non-zero, output the new number of kpoints on unit iout
!! kbz(3,nkbz)= k vectors in the BZ.
!! nkbz = number of k-points whose weights are wtk
!! nsym=number of space group symmetries
!! symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!! timrev: if 1, the time reversal operation has to be taken into account
!!         if 0, no time reversal symmetry.
!! wtk(nkbz)=weight assigned to each k point.
!! comm=MPI communicator.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine prt_eigocc(eigen,eltemp,iout,mband,nband,nkpt,nsppol,occ,prefix_out,ucvol,wtk)

!Arguments -------------------------------
!scalars
  integer,intent(in) :: iout,mband,nkpt,nsppol
  real(dp),intent(in) :: eltemp,ucvol
  character(len=*),intent(in) :: prefix_out
!arrays
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: wtk(nkpt)


!Local variables -------------------------
!scalars
  integer :: ii,ikpt,isppol,nband_k,temp_unit
  character(len=200) :: fname_eigocc
  character(len=500) :: msg

! *********************************************************************

  fname_eigocc = trim(prefix_out) // trim('_EIGOCC')
  write(msg,'(a,a)') ' prt_eigocc : about to open file ',trim(fname_eigocc)
  call wrtout(iout,msg,'COLL')
  if (open_file(fname_eigocc,msg,newunit=temp_unit,status='unknown',form='formatted') /= 0) then
    MSG_ERROR(msg)
  end if
  rewind(temp_unit)

! Loop over spins
  do isppol=1,nsppol 
  
!   Loop over k-points
    do ikpt=1,nkpt
      nband_k=nband(ikpt+(isppol-1)*nkpt)
      
!     Loop over bands
      do ii=0,(nband_k-1)/4

      end do
    end do
  end do
  

end subroutine prt_eigocc

end module hightemp
