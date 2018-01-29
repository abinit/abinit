!{\src2tex{textfont=tt}}
!!****f* ABINIT/getghcnd
!!
!! NAME
!! getghcnd
!!
!! FUNCTION
!! Compute <G|H_ND|C> for input vector |C> expressed in reciprocal space
!! Result is put in array ghcnc. H_ND is the Hamiltonian due to magnetic dipoles
!! on the nuclear sites.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cwavef(2,npw*nspinor*ndat)=planewave coefficients of wavefunction.
!! gs_ham <type(gs_hamiltonian_type)>=all data for the Hamiltonian to be applied
!! my_nspinor=number of spinorial components of the wavefunctions (on current proc)
!! ndat=number of FFT to do in //
!!
!! OUTPUT
!! ghcnd(2,npw*my_nspinor*ndat)=matrix elements <G|H_ND|C>
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Application of <k^prime|H|k> or <k|H|k^prime> not implemented!
!!
!! PARENTS
!!      getghc
!!
!! CHILDREN
!!      zhpmv
!!
!! NOTES
!!  This routine applies the Hamiltonian due to an array of magnetic dipoles located
!!  at the atomic nuclei to the input wavefunction. Strategy below is to take advantage of
!!  Hermiticity to store H_ND in triangular form and then use a BLAS call to zhpmv to apply to
!!  input vector in one shot.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine getghcnd(cwavef,ghcnd,gs_ham,my_nspinor,ndat)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

 use m_hamiltonian, only : gs_hamiltonian_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getghcnd'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_nspinor,ndat
 type(gs_hamiltonian_type),intent(in),target :: gs_ham
!arrays
 real(dp),intent(in) :: cwavef(2,gs_ham%npw_k*my_nspinor*ndat)
 real(dp),intent(out) :: ghcnd(2,gs_ham%npw_k*my_nspinor*ndat)

!Local variables-------------------------------
!scalars
 integer :: cwavedim,igp
 character(len=500) :: message
 !arrays
 complex(dpc),allocatable :: inwave(:),hggc(:)

! *********************************************************************

 if (gs_ham%matblk /= gs_ham%natom) then
   write(message,'(a,i4,a,i4)')' gs_ham%matblk = ',gs_ham%matblk,' but natom = ',gs_ham%natom
   MSG_ERROR(message)
 end if
 if (ndat /= 1) then
   write(message,'(a,i4,a)')' ndat = ',ndat,' but getghcnd requires ndat = 1'
   MSG_ERROR(message)
 end if
 if (my_nspinor /= 1) then
   write(message,'(a,i4,a)')' nspinor = ',my_nspinor,' but getghcnd requires nspinor = 1'
   MSG_ERROR(message)
 end if
 if (any(abs(gs_ham%kpt_k(:)-gs_ham%kpt_kp(:))>tol8)) then
   message=' not allowed for kpt(left)/=kpt(right)!'
   MSG_BUG(message)
 end if

 cwavedim = gs_ham%npw_k*my_nspinor*ndat
 ABI_ALLOCATE(hggc,(cwavedim))
 ABI_ALLOCATE(inwave,(cwavedim))

 do igp = 1, gs_ham%npw_k
   inwave(igp) = cmplx(cwavef(1,igp),cwavef(2,igp),kind=dpc)
 end do

! apply hamiltonian hgg to input wavefunction inwave, result in hggc
 call ZHPMV('L',cwavedim,cone,gs_ham%nucdipmom_k,inwave,1,czero,hggc,1)

 do igp=1,gs_ham%npw_k
   ghcnd(1,igp) = dreal(hggc(igp))
   ghcnd(2,igp) = dimag(hggc(igp))
 end do

 ABI_DEALLOCATE(hggc)
 ABI_DEALLOCATE(inwave)

! ghcnd=zero

 ! ABI_DEALLOCATE(hgg)

end subroutine getghcnd
!!***
