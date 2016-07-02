!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffReadSkipK
!! NAME
!! WffReadSkipK
!!
!! FUNCTION
!!  (Wavefunction file, read action : skip one k-point blok)
!!  This subroutine skips the block of records
!!  related to one k point, and one spin-polarization, that
!!  contains the wavefunctions as well as the eigenvalues and occupations,
!!  in a wavefunction file that has been already initialized.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  formeig=format of the eigenvalues
!!   0 => vector of eigenvalues (for Ground-State files)
!!   1 => hermitian matrix of eigenvalues (for Response-Function files)
!!  headform=format of the header of the wf file, also governing the k block format
!!   in case headform=0, use the default (current) format and headform
!!  ikpt=index of current k point (only needed for error message)
!!  isppol=spin polarization currently treated (only needed for error message)
!!  mpi_enreg=informations about MPI parallelization
!!  wff=structured info for wavefunction file
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      d2frnl,dfpt_nstdy,dfpt_nstpaw,dfpt_vtorho,initwf,newkpt,wfsinp
!!
!! CHILDREN
!!      rwwf
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine WffReadSkipK(formeig,headform,ikpt,isppol,mpi_enreg,wff)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_wffile

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'WffReadSkipK'
 use interfaces_56_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: formeig,headform,ikpt,isppol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
!scalars
 integer :: icg,mband,mcg,nband,nband_disk,option,optkg,tim_rwwf
 integer,parameter :: nspinor1=1,npw1=1
!arrays
 integer,allocatable :: kg_dum(:,:)
 real(dp) :: cg_dum(2,1),occ_dum(1)
 real(dp),allocatable :: eig_dum(:)

! *************************************************************************

 ! No need to skip if netcdf
 if (wff%iomode == IO_MODE_ETSF) return

 option=-1
 tim_rwwf=0 ; mcg=1 ; mband=1 ; icg=0 ; optkg=0 ; nband=0
 ABI_ALLOCATE(eig_dum,(2**formeig))
 ABI_ALLOCATE(kg_dum,(3,optkg*npw1))

 call rwwf(cg_dum,eig_dum,formeig,headform,icg,ikpt,isppol,kg_dum,mband,mcg,mpi_enreg,nband,&
& nband_disk,npw1,nspinor1,occ_dum,option,optkg,tim_rwwf,wff)

 ABI_DEALLOCATE(eig_dum)
 ABI_DEALLOCATE(kg_dum)

end subroutine WffReadSkipK
!!***
