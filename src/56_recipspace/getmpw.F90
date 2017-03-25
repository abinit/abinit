!{\src2tex{textfont=tt}}
!!****f* ABINIT/getmpw
!! NAME
!! getmpw
!!
!! FUNCTION
!! From input ecut, combined with ucvol and gmet, compute recommended mpw
!! mpw is the maximum number of plane-waves in the wave-function basis
!!  for one processor of the WF group
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! ecut=plane-wave cutoff energy in Hartrees
!! exchn2n3d=if n1, n2 and n3 are exchanged
!! gmet(3,3)=reciprocal space metric (bohr**-2).
!! istwfk(nkpt)=input option parameter that describes the storage of wfs
!! kptns(3,nkpt)=real(dp) array for k points (normalisation is already taken into account)
!! mpi_enreg=information about MPI parallelization
!! nkpt=integer number of k points in the calculation
!!
!! OUTPUT
!! mpw=maximal number of plane waves over all k points of the processor
!!  (for one processor of the WF group)
!!
!! PARENTS
!!      dfpt_looppert,memory_eval,mpi_setup,scfcv
!!
!! CHILDREN
!!      kpgsph,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine getmpw(ecut,exchn2n3d,gmet,istwfk,kptns,mpi_enreg,mpw,nkpt)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

 use m_fftcore, only : kpgsph

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getmpw'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: exchn2n3d,nkpt
 integer,intent(out) :: mpw
 real(dp),intent(in) :: ecut
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: istwfk(nkpt)
 real(dp),intent(in) :: gmet(3,3),kptns(3,nkpt)

!Local variables-------------------------------
!scalars
 integer :: ikpt,istwf_k,npw
! integer :: npwwrk,pad=50
! real(dp) :: scale=1.3_dp
 character(len=500) :: message
!arrays
 integer,allocatable :: kg(:,:)
 real(dp) :: kpoint(3)

! *************************************************************************

!An upper bound for mpw, might be obtained as follows
!the average number of plane-waves in the cutoff sphere is
!npwave = (2*ecut)**(3/2) * ucvol / (6*pi**2)
!the upper bound is calculated as
!npwwrk = int(scale * npwave) + pad
!rescale so an upper bound
!npwave=nint(ucvol*(2.0_dp*ecut)**1.5_dp/(6.0_dp*pi**2))
!npwwrk=nint(dble(npwave)*scale)+pad

 ABI_ALLOCATE(kg,(3,100))

!set mpw to zero, as needed for only counting in kpgsph
 mpw = 0

!Might be parallelized over k points ? !
 do ikpt = 1,nkpt
!  Do computation of G sphere, returning npw
   kpoint(:)=kptns(:,ikpt)
   istwf_k=istwfk(ikpt)
   call kpgsph(ecut,exchn2n3d,gmet,0,ikpt,istwf_k,kg,kpoint,0,mpi_enreg,0,npw)
   mpw = max(npw,mpw)
 end do

 write(message,'(a,i0)') ' getmpw: optimal value of mpw= ',mpw
 call wrtout(std_out,message,'COLL')

 ABI_DEALLOCATE(kg)

end subroutine getmpw
!!***
