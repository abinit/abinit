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
!! Copyright (C) 1998-2017 ABINIT group
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
 integer :: row,col,coldx,ndp_index,npw
 real(dp) :: Aijr, Aiji,coli,colr,rowi,rowr
 character(len=500) :: message
 !arrays

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

! size of incoming wavefunction
 ! cwavedim = gs_ham%npw_k*my_nspinor*ndat
 ! ABI_ALLOCATE(hgg,(cwavedim*(cwavedim+1)/2))
 ! ABI_ALLOCATE(hggc,(cwavedim))
 ! ABI_ALLOCATE(inwave,(cwavedim))

 ! Hand-coded H_ND.cwave (715 sec) Least memory, ok performance
 ghcnd=zero
 npw = gs_ham%npw_k
 ndp_index = 0
 do col=1, npw
    colr = cwavef(1,col)
    coli = cwavef(2,col)
    coldx = (2*npw-col+2)*(col-1)/2 ! elements of columns already completed, each
    ! one less planewave than the one before, starting with npw in col 1
    do row = col+1, npw  ! row == col term is strictly zero
       ndp_index = coldx + row - col + 1
       Aijr = gs_ham%nucdipmom_k(1,ndp_index)
       Aiji = gs_ham%nucdipmom_k(2,ndp_index)
       rowr = cwavef(1,row)
       rowi = cwavef(2,row)
       ghcnd(1,row) = ghcnd(1,row) + Aijr*colr - Aiji*coli
       ghcnd(2,row) = ghcnd(2,row) + Aijr*coli + Aiji*colr
       ghcnd(1,col) = ghcnd(1,col) + Aijr*rowr + Aiji*rowi
       ghcnd(2,col) = ghcnd(2,col) + Aijr*rowi - Aiji*rowr
    end do
 end do

 ! MATMUL with real objects (1631 sec) most memory, worst performance
 ! cwavedim = gs_ham%npw_k*my_nspinor*ndat
 ! ghcnd(1:2,1:cwavedim)=zero
 ! ABI_ALLOCATE(hreal,(cwavedim,cwavedim))
 ! ABI_ALLOCATE(himag,(cwavedim,cwavedim))
 ! ABI_ALLOCATE(cwavereal,(cwavedim))
 ! ABI_ALLOCATE(cwaveimag,(cwavedim))
 ! cwavereal(1:cwavedim) = cwavef(1,1:cwavedim)
 ! cwaveimag(1:cwavedim) = cwavef(2,1:cwavedim)
 ! ndp_index = 0
 ! do col=1, gs_ham%npw_k
 !    do row = col, gs_ham%npw_k 
 !       ndp_index = ndp_index + 1
 !       hreal(row,col) = gs_ham%nucdipmom_k(1,ndp_index)
 !       hreal(col,row) = hreal(row,col)
 !       himag(row,col) = gs_ham%nucdipmom_k(2,ndp_index)
 !       himag(col,row) = -himag(row,col)
 !    end do
 ! end do
 ! ghcnd(1,1:cwavedim) = MATMUL(hreal,cwavereal) - MATMUL(himag,cwaveimag)
 ! ghcnd(2,1:cwavedim) = MATMUL(hreal,cwaveimag) + MATMUL(himag,cwavereal)
 ! ABI_DEALLOCATE(hreal)
 ! ABI_DEALLOCATE(himag)
 ! ABI_DEALLOCATE(cwavereal)
 ! ABI_DEALLOCATE(cwaveimag)

 ! BLAS with real objects (1401 sec) most memory, poor performance
 ! cwavedim = gs_ham%npw_k*my_nspinor*ndat
 ! ghcnd(1:2,1:cwavedim)=zero
 ! ABI_ALLOCATE(hreal,(cwavedim,cwavedim))
 ! ABI_ALLOCATE(himag,(cwavedim,cwavedim))
 ! ABI_ALLOCATE(cwavereal,(cwavedim))
 ! ABI_ALLOCATE(cwaveimag,(cwavedim))
 ! cwavereal(1:cwavedim) = cwavef(1,1:cwavedim)
 ! cwaveimag(1:cwavedim) = cwavef(2,1:cwavedim)
 ! ndp_index = 0
 ! do col=1, gs_ham%npw_k
 !    do row = col, gs_ham%npw_k 
 !       ndp_index = ndp_index + 1
 !       hreal(row,col) = gs_ham%nucdipmom_k(1,ndp_index)
 !       hreal(col,row) = hreal(row,col)
 !       himag(row,col) = gs_ham%nucdipmom_k(2,ndp_index)
 !       himag(col,row) = -himag(row,col)
 !    end do
 ! end do

 ! ABI_ALLOCATE(work,(cwavedim))
 ! call DGEMV('N',cwavedim,cwavedim,1.0D0,hreal,cwavedim,cwavereal,1,0.0D0,work,1)
 ! call DGEMV('N',cwavedim,cwavedim,-1.0D0,himag,cwavedim,cwaveimag,1,1.0D0,work,1)
 ! ghcnd(1,1:cwavedim) = work(1:cwavedim)

 ! call DGEMV('N',cwavedim,cwavedim,1.0D0,hreal,cwavedim,cwaveimag,1,0.0D0,work,1)
 ! call DGEMV('N',cwavedim,cwavedim,1.0D0,himag,cwavedim,cwavereal,1,1.0D0,work,1)
 ! ghcnd(2,1:cwavedim) = work(1:cwavedim)
 
 ! ABI_DEALLOCATE(hreal)
 ! ABI_DEALLOCATE(himag)
 ! ABI_DEALLOCATE(cwavereal)
 ! ABI_DEALLOCATE(cwaveimag)
 ! ABI_DEALLOCATE(work)

 ! MATMUL with complex objects (1353 sec) much memory, poor performance
 ! cwavedim = gs_ham%npw_k*my_nspinor*ndat
 ! ABI_ALLOCATE(hgg,(cwavedim,cwavedim))
 ! ABI_ALLOCATE(inwave,(cwavedim))
 ! ABI_ALLOCATE(work,(cwavedim))
 ! ndp_index = 0
 ! do col=1, gs_ham%npw_k
 !    inwave(col)=CMPLX(cwavef(1,col),cwavef(2,col),kind=dpc)
 !    do row = col, gs_ham%npw_k 
 !       ndp_index = ndp_index + 1
 !       hgg(row,col) = CMPLX(gs_ham%nucdipmom_k(1,ndp_index),gs_ham%nucdipmom_k(2,ndp_index),kind=dpc)
 !       hgg(col,row) = CONJG(hgg(row,col))
 !    end do
 ! end do

 ! work = MATMUL(hgg,inwave)
 ! do col = 1, cwavedim
 !    ghcnd(1,col) = REAL(work(col))
 !    ghcnd(2,col) = AIMAG(work(col))
 ! end do
 

 ! ABI_DEALLOCATE(hgg)
 ! ABI_DEALLOCATE(inwave)
 ! ABI_DEALLOCATE(work)

 ! BLAS with complex objects (656 sec) much memory, best performance
 ! cwavedim = gs_ham%npw_k*my_nspinor*ndat
 ! cwavedimp = cwavedim*(cwavedim+1)/2
 ! ABI_ALLOCATE(hgg,(cwavedimp))
 ! ABI_ALLOCATE(inwave,(cwavedim))
 ! ABI_ALLOCATE(work,(cwavedim))
 ! ndp_index = 0
 ! do col=1, gs_ham%npw_k
 !    inwave(col)=CMPLX(cwavef(1,col),cwavef(2,col),kind=dpc)
 !    do row = col, gs_ham%npw_k 
 !       ndp_index = ndp_index + 1
 !       hgg(ndp_index) = CMPLX(gs_ham%nucdipmom_k(1,ndp_index),gs_ham%nucdipmom_k(2,ndp_index),kind=dpc)
 !    end do
 ! end do

 ! call zhpmv('L',cwavedim,cone,hgg,inwave,1,czero,work,1)

 ! do col = 1, cwavedim
 !    ghcnd(1,col) = REAL(work(col))
 !    ghcnd(2,col) = AIMAG(work(col))
 ! end do
 
 
 ! ABI_DEALLOCATE(hgg)
 ! ABI_DEALLOCATE(inwave)
 ! ABI_DEALLOCATE(work)

!  do igp = 1, gs_ham%npw_k
!    inwave(igp) = cmplx(cwavef(1,igp),cwavef(2,igp),kind=dpc)
!  end do

!  do igp = 1, gs_ham%npw_k*(gs_ham%npw_k+1)/2
!     hgg(igp) = cmplx(gs_ham%nucdipmom_k(1,igp),gs_ham%nucdipmom_k(2,igp),kind=dpc)
!  end do

! ! apply hamiltonian hgg to input wavefunction inwave, result in hggc
!  call zhpmv('L',cwavedim,cone,hgg,inwave,1,czero,hggc,1)

!  do igp=1,gs_ham%npw_k
!    ghcnd(1,igp) = dreal(hggc(igp))
!    ghcnd(2,igp) = dimag(hggc(igp))
!  end do

! ghcnd=zero

 ! ABI_DEALLOCATE(hgg)
 ! ABI_DEALLOCATE(hggc)
 ! ABI_DEALLOCATE(inwave)

end subroutine getghcnd
!!***
