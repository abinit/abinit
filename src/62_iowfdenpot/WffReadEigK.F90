!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffReadEigK
!! NAME
!! WffReadEigK
!!
!! FUNCTION
!!  (Wavefunction file, read eigenvalues at one k point)
!!  This subroutine reads the block of records
!!  related to one k point, and one spin-polarization, that
!!  contains the wavefunctions as well as the eigenvalues and occupations .
!!  It outputs only the eigenvalues.
!!  The wavefunction file should have been prepared
!!  outside of this routine, in order to read or write the correct records.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  formeig=format of the eigenvalues
!!   0 => vector of eigenvalues (Ground-State case)
!!   1 => hermitian matrix of eigenvalues (Response Function case)
!!  headform=format of the header of the wf file, also governing the k block format
!!   in case headform=0, use the default (current) format and headform
!!  ikpt=index of current k point (only needed for error message)
!!  isppol=spin polarization currently treated (only needed for error message)
!!  mband=maximum number of bands (governs the dimension of eigen)
!!  mpi_enreg=information about MPI parallelization
!!  nband=number of bands actually in eigen
!!    can be equal, larger or smaller than nband_disk, but
!!    eigen will not be completely filled if nband>nband_disk)
!!    must be less or equal to mband
!!  tim_rwwf=timing code of the calling routine (set to 0 if not attributed)
!!  wff=structured info for the wavefunction file
!!
!! OUTPUT
!!  eigen((2*mband)**formeig *mband)=array for holding eigenvalues (hartree)
!!
!! PARENTS
!!
!! CHILDREN
!!      rwwf
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine WffReadEigK(eigen,formeig,headform,ikpt,isppol,mband,mpi_enreg,nband,tim_rwwf,wff)

 use defs_basis
 use defs_abitypes
 use m_wffile
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'WffReadEigK'
 use interfaces_56_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: formeig,headform,ikpt,isppol,mband,nband,tim_rwwf
 type(MPI_type),intent(in) :: mpi_enreg
 type(wffile_type),intent(inout) :: wff
!arrays
 real(dp),intent(out) :: eigen((2*mband)**formeig*mband)

!Local variables-------------------------------
!scalars
 integer,parameter :: nspinor1=1
 integer :: icg,mcg,nband_disk,npw,option,optkg
!arrays
 integer,allocatable :: kg_dum(:,:)
 real(dp) :: cg_dum(2,1)
 real(dp),allocatable :: occ_dum(:)

! *************************************************************************

!DEBUG
!write(std_out,*)' WffReadEigK : enter '
!ENDDEBUG

 option=3
 mcg=1 ; icg=0 ; optkg=0 
 nband_disk=mband; npw=0
 ABI_ALLOCATE(kg_dum,(3,optkg))
 ABI_ALLOCATE(occ_dum,(mband))

!DEBUG
!write(std_out,*)' WffReadEigK : formeig,option=',formeig,option
!ENDDEBUG

 call rwwf(cg_dum,eigen,formeig,headform,icg,ikpt,isppol,kg_dum,mband,mcg,mpi_enreg,&
& nband,nband_disk,&
& npw,nspinor1,occ_dum,option,optkg,tim_rwwf,wff)

 ABI_DEALLOCATE(kg_dum)
 ABI_DEALLOCATE(occ_dum)

end subroutine WffReadEigK
!!***
