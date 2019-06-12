
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"
!
! Copyright (C) 2014 Mitsuaki Kawamura
!
! Permission is hereby granted, free of charge, to any person obtaining a
! copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
! 
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
! TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
MODULE libtetrabz_occ_mod
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: libtetrabz_occ, libtetrabz_fermieng
  !
CONTAINS
!
! Compute occupation
!
SUBROUTINE libtetrabz_occ(ltetra,bvec,nb,nge,eig,ngw,wght,comm) 
  !
  USE ISO_C_BINDING
  USE libtetrabz_common, ONLY : libtetrabz_initialize, libtetrabz_interpol_indx, libtetrabz_mpisum_dv
  IMPLICIT NONE
  !
  INTEGER(C_INT),INTENT(IN) :: ltetra, nb, nge(3), ngw(3)
  REAL(C_DOUBLE),INTENT(IN) :: bvec(9), eig(nb,PRODUCT(nge(1:3)))
  REAL(C_DOUBLE),INTENT(OUT) :: wght(nb,PRODUCT(ngw(1:3)))
  INTEGER(C_INT),INTENT(IN),OPTIONAL :: comm
  !
  LOGICAL :: linterpol
  INTEGER :: nt_local, nk_local, nkBZ, ik, kintp(20), nintp
  INTEGER,ALLOCATABLE :: ik_global(:,:), ik_local(:,:)
  REAL(8) :: wlsm(4,20), wintp(1,20)
  REAL(8),ALLOCATABLE :: wghtd(:,:,:), kvec(:,:)
  !
  nintp = 16 * ltetra - 12
  !
  IF(PRESENT(comm)) THEN
     CALL libtetrabz_initialize(ltetra,nge,ngw,bvec,linterpol,wlsm,nk_local,&
     &                          nt_local,nkBZ,ik_global,ik_local,kvec,comm)
  ELSE
     CALL libtetrabz_initialize(ltetra,nge,ngw,bvec,linterpol,wlsm,nk_local,&
     &                          nt_local,nkBZ,ik_global,ik_local,kvec)
  END IF
  !
  IF(linterpol) THEN
     !
     ABI_MALLOC(wghtd, (nb,1,nk_local))
     CALL libtetrabz_occ_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig,nk_local,wghtd,0d0)
     !
     ! Interpolation
     !
     wght(1:nb,1:PRODUCT(ngw(1:3))) = 0d0
     DO ik = 1, nk_local
        CALL libtetrabz_interpol_indx(nintp,ngw,kvec(1:3,ik),kintp,wintp)
        wght(1:nb,kintp(1:nintp)) = wght(1:nb,             kintp(1:nintp)) &
        &                 + MATMUL(wghtd(1:nb,1:1,ik), wintp(1:1,1:nintp))
     END DO ! ik = 1, nk_local
     ABI_FREE(wghtd)
     ABI_FREE(kvec)
     !
     IF(PRESENT(comm)) CALL libtetrabz_mpisum_dv(comm, nb * PRODUCT(ngw(1:3)), wght)
     !
  ELSE
     CALL libtetrabz_occ_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig,nk_local,wght,0d0)
  END IF
  !
  ABI_FREE(ik_global)
  ABI_FREE(ik_local)
  !
END SUBROUTINE libtetrabz_occ
!
! Calculate Fermi energy
!
SUBROUTINE libtetrabz_fermieng(ltetra,bvec,nb,nge,eig,ngw,wght,ef,nelec,comm) 
  !
  USE ISO_C_BINDING
  USE libtetrabz_common, ONLY : libtetrabz_initialize, libtetrabz_interpol_indx, &
  &                             libtetrabz_mpisum_d, libtetrabz_mpisum_dv
  IMPLICIT NONE
  !
  INTEGER(C_INT),INTENT(IN) :: ltetra, nb, nge(3), ngw(3)
  REAL(C_DOUBLE),INTENT(IN) :: bvec(9), nelec, eig(nb,PRODUCT(nge(1:3)))
  REAL(C_DOUBLE),INTENT(OUT) :: ef
  REAL(C_DOUBLE),INTENT(OUT) :: wght(nb,PRODUCT(ngw(1:3)))
  INTEGER(C_INT),INTENT(IN),OPTIONAL :: comm
  !
  LOGICAL :: linterpol
  INTEGER :: nt_local, nk_local, nkBZ, ik, kintp(20), nintp
  INTEGER,ALLOCATABLE :: ik_global(:,:), ik_local(:,:)
  REAL(8) :: wlsm(4,20), wintp(1,20)
  REAL(8),ALLOCATABLE :: wghtd(:,:,:), kvec(:,:)
  !
  INTEGER :: iter, maxiter = 300
  REAL(8) :: elw, eup, eps= 1d-10, sumkmid
  !
  nintp = 16 * ltetra - 12
  !
  IF(PRESENT(comm)) THEN
     CALL libtetrabz_initialize(ltetra,nge,ngw,bvec,linterpol,wlsm,nk_local,&
     &                          nt_local,nkBZ,ik_global,ik_local,kvec,comm)
  ELSE
     CALL libtetrabz_initialize(ltetra,nge,ngw,bvec,linterpol,wlsm,nk_local,&
     &                          nt_local,nkBZ,ik_global,ik_local,kvec)
  END IF
  !
  IF(linterpol) then
    ABI_MALLOC(wghtd, (nb,1,nk_local))
  end if
  !
  elw = MINVAL(eig(1:nb,1:PRODUCT(nge(1:3))))
  eup = MAXVAL(eig(1:nb,1:PRODUCT(nge(1:3))))
  !
  ! Bisection method
  !
  DO iter = 1, maxiter
     !
     ef = (eup + elw) / 2.d0
     !
     ! Calc. # of electrons 
     !
     IF(linterpol) THEN
        CALL libtetrabz_occ_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig,nk_local,wghtd,ef)
        sumkmid = SUM(wghtd(1:nb,1,1:nk_local))
        !
        IF(PRESENT(comm)) CALL libtetrabz_mpisum_d(comm,sumkmid)
        !
     ELSE
        CALL libtetrabz_occ_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig,nk_local,wght,ef)
        sumkmid = SUM(wght(1:nb,1:nk_local))
     END IF
     !
     ! convergence check
     !
     IF(ABS(sumkmid - nelec) < eps) THEN
        EXIT
     ELSE IF(sumkmid < nelec) THEN
        elw = ef
     ELSE
        eup = ef
     ENDIF
     !
  END DO  ! iter
  !
  IF(iter >= maxiter) STOP "libtetrabz_fermieng"
  !
  IF(linterpol) THEN
     !
     ! Interpolation
     !
     wght(1:nb,1:PRODUCT(ngw(1:3))) = 0d0
     DO ik = 1, nk_local
        CALL libtetrabz_interpol_indx(nintp,ngw,kvec(1:3,ik),kintp,wintp)
        wght(1:nb,kintp(1:nintp)) = wght(1:nb,             kintp(1:nintp)) &
        &                 + MATMUL(wghtd(1:nb,1:1,ik), wintp(1:1,1:nintp))
     END DO ! ik = 1, nk_local
     ABI_FREE(wghtd)
     ABI_FREE(kvec)
     !
     IF(PRESENT(comm)) CALL libtetrabz_mpisum_dv(comm, nb * PRODUCT(ngw(1:3)), wght)
     !
  END IF ! (linterpol)
  !
  ABI_FREE(ik_global)
  ABI_FREE(ik_local)
  !
END SUBROUTINE libtetrabz_fermieng
!
! Main SUBROUTINE for occupation : Theta(EF - E1)
!
SUBROUTINE libtetrabz_occ_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig,nk_local,occ,ef)
  !
  USE libtetrabz_common, ONLY : libtetrabz_sort, &
  &                             libtetrabz_tsmall_a1, libtetrabz_tsmall_b1, &
  &                             libtetrabz_tsmall_b2, libtetrabz_tsmall_b3, &
  &                             libtetrabz_tsmall_c1, libtetrabz_tsmall_c2, &
  &                             libtetrabz_tsmall_c3
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: nt_local, nb, nkBZ, nk_local, &
  &                     ik_global(20,nt_local), ik_local(20,nt_local)
  REAL(8),INTENT(IN) :: wlsm(4,20), ef, eig(nb,nkBZ)
  REAL(8),INTENT(OUT) :: occ(nb,nk_local)
  !
  INTEGER :: ib, indx(4), it
  REAL(8) :: e(4), ei1(4,nb), tsmall(4,4), V, w1(4)
  !
  occ(1:nb,1:nk_local) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(ef,eig,ik_global,ik_local,nb,nt_local,occ,wlsm) &
  !$OMP & PRIVATE(e,ei1,ib,indx,it,tsmall,V,w1)
  !
  DO it = 1, nt_local
     !
     DO ib = 1, nb
        ei1(1:4,ib) = MATMUL(wlsm(1:4,1:20), eig(ib,ik_global(1:20,it)))
     END DO
     !
     !$OMP DO
     DO ib = 1, nb
        !
        w1(1:4) = 0d0
        e(1:4) = ei1(1:4, ib) - ef
        CALL libtetrabz_sort(4,e,indx)
        !
        IF(e(1) <= 0d0 .AND. 0d0 < e(2)) THEN
           !
           CALL libtetrabz_tsmall_a1(e,V,tsmall)
           w1(indx(1:4)) = w1(indx(1:4)) + V * SUM(tsmall(1:4,1:4), 1) * 0.25d0
           !
        ELSE IF(e(2) <= 0d0 .AND. 0d0 < e(3)) THEN
           !
           CALL libtetrabz_tsmall_b1(e,V,tsmall)
           w1(indx(1:4)) = w1(indx(1:4)) + V * SUM(tsmall(1:4,1:4), 1) * 0.25d0
           !
           CALL libtetrabz_tsmall_b2(e,V,tsmall)
           w1(indx(1:4)) = w1(indx(1:4)) + V * SUM(tsmall(1:4,1:4), 1) * 0.25d0
           !
           CALL libtetrabz_tsmall_b3(e,V,tsmall)
           w1(indx(1:4)) = w1(indx(1:4)) + V * SUM(tsmall(1:4,1:4), 1) * 0.25d0
           !
        ELSE IF(e(3) <= 0d0 .AND. 0d0 < e(4)) THEN
           !
           CALL libtetrabz_tsmall_c1(e,V,tsmall)
           w1(indx(1:4)) = w1(indx(1:4)) + V * SUM(tsmall(1:4,1:4), 1) * 0.25d0
           !
           CALL libtetrabz_tsmall_c2(e,V,tsmall)
           w1(indx(1:4)) = w1(indx(1:4)) + V * SUM(tsmall(1:4,1:4), 1) * 0.25d0
           !
           CALL libtetrabz_tsmall_c3(e,V,tsmall)
           w1(indx(1:4)) = w1(indx(1:4)) + V * SUM(tsmall(1:4,1:4), 1) * 0.25d0
           !
        ELSE IF(e(4) <= 0d0) THEN
           !
           w1(1:4) = 0.25d0
           !
        ELSE
           !
           CYCLE
           !
        END IF
        !
        occ(ib,ik_local(1:20,it)) = occ(ib,ik_local(1:20,it)) &
        &                + MATMUL(w1(1:4), wlsm(1:4,1:20))
        !
     END DO ! ib
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !
  !$OMP END PARALLEL
  !
  occ(1:nb,1:nk_local) = occ(1:nb,1:nk_local) / DBLE(6 * nkBZ)
  !
END SUBROUTINE libtetrabz_occ_main
!
END MODULE libtetrabz_occ_mod
