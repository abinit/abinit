
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
MODULE libtetrabz_dbldelta_mod
  !

  use m_abicore
  use m_errors

  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: libtetrabz_dbldelta
  !
CONTAINS
!
! Compute doubledelta
!
SUBROUTINE libtetrabz_dbldelta(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,comm)
  !
  USE ISO_C_BINDING
  USE libtetrabz_common, ONLY : libtetrabz_initialize, libtetrabz_interpol_indx, libtetrabz_mpisum_dv
  IMPLICIT NONE
  !
  INTEGER(C_INT),INTENT(IN) :: ltetra, nb, nge(3), ngw(3)
  REAL(C_DOUBLE),INTENT(IN) :: bvec(9), eig1(nb,PRODUCT(nge(1:3))), eig2(nb,PRODUCT(nge(1:3)))
  REAL(C_DOUBLE),INTENT(OUT) :: wght(nb*nb,PRODUCT(ngw(1:3)))
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
     ABI_MALLOC(wghtd, (nb*nb,1,nk_local))
     CALL libtetrabz_dbldelta_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig1,eig2,nk_local,wghtd)
     !
     ! Interpolation
     !
     wght(1:nb*nb,1:PRODUCT(ngw(1:3))) = 0d0
     DO ik = 1, nk_local
        CALL libtetrabz_interpol_indx(nintp,ngw,kvec(1:3,ik),kintp,wintp)
        wght(1:nb*nb,kintp(1:nintp)) = wght(1:nb*nb,             kintp(1:nintp)) &
        &                    + MATMUL(wghtd(1:nb*nb,1:1,ik), wintp(1:1,1:nintp))
     END DO ! ik = 1, nk_local
     ABI_FREE(wghtd)
     ABI_FREE(kvec)
     !
     IF(PRESENT(comm)) CALL libtetrabz_mpisum_dv(comm, nb*nb*PRODUCT(ngw(1:3)), wght)
     !
  ELSE
     CALL libtetrabz_dbldelta_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig1,eig2,nk_local,wght)
  END IF
  !
  ABI_FREE(ik_global)
  ABI_FREE(ik_local)
  !
END SUBROUTINE libtetrabz_dbldelta
!
! Main SUBROUTINE for Delta(E1) * Delta(E2)
!
SUBROUTINE libtetrabz_dbldelta_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig1,eig2,nk_local,dbldelta)
  !
  USE libtetrabz_common, ONLY : libtetrabz_sort, &
  &                             libtetrabz_triangle_a1, libtetrabz_triangle_b1, &
  &                             libtetrabz_triangle_b2, libtetrabz_triangle_c1
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: nt_local, nb, nkBZ, nk_local, &
  &                     ik_global(20,nt_local), ik_local(20,nt_local)
  REAL(8),INTENT(IN) :: wlsm(4,20), eig1(nb,nkBZ), eig2(nb,nkBZ)
  REAL(8),INTENT(OUT) :: dbldelta(nb,nb,nk_local)
  !
  INTEGER :: ib, indx(4), it
  REAL(8) :: e(4), ei1(4,nb), ej1(4,nb), ej2(3,nb), V, thr = 1d-10, &
  &          tsmall(3,4), w1(nb,4), w2(nb,3)
  !
  dbldelta(1:nb,1:nb,1:nk_local) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(dbldelta,eig1,eig2,ik_global,ik_local,nb,nt_local,thr,wlsm) &
  !$OMP & PRIVATE(e,ei1,ej1,ej2,ib,indx,it,tsmall,V,w1,w2)
  !
  DO it = 1, nt_local
     !
     DO ib = 1, nb
        ei1(1:4,ib) = MATMUL(wlsm(1:4,1:20), eig1(ib, ik_global(1:20,it)))
        ej1(1:4,ib) = MATMUL(wlsm(1:4,1:20), eig2(ib, ik_global(1:20,it)))
     END DO
     !
     !$OMP DO
     DO ib = 1, nb
        !
        w1(1:nb,1:4) = 0d0
        e(1:4) = ei1(1:4, ib)
        CALL libtetrabz_sort(4,e,indx)
        !
        IF(e(1) < 0d0 .AND. 0d0 <= e(2)) THEN
           !
           CALL libtetrabz_triangle_a1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ej2(1:3,1:nb) = MATMUL(tsmall(1:3,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_dbldelta2(nb,ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:3), tsmall(1:3,1:4))
              !
           END IF
           !
        ELSE IF( e(2) < 0d0 .AND. 0d0 <= e(3)) THEN
           !
           CALL libtetrabz_triangle_b1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ej2(1:3,1:nb) = MATMUL(tsmall(1:3,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_dbldelta2(nb,ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:3), tsmall(1:3,1:4))
              !
           END IF
           !
           CALL libtetrabz_triangle_b2(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ej2(1:3,1:nb) = MATMUL(tsmall(1:3,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_dbldelta2(nb,ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:3), tsmall(1:3,1:4))
              !
           END IF
           !
        ELSE IF(e(3) < 0d0 .AND. 0d0 < e(4)) THEN
           !
           CALL libtetrabz_triangle_c1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ej2(1:3,1:nb) = MATMUL(tsmall(1:3,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_dbldelta2(nb,ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:3), tsmall(1:3,1:4))
              !
           END IF
           !
        END IF
        !
        dbldelta(1:nb,ib,ik_local(1:20,it)) = dbldelta(1:nb,ib,   ik_local(1:20,it)) &
        &                                  + MATMUL(w1(1:nb,1:4), wlsm(1:4,1:20))
        !
     END DO ! ib
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !
  !$OMP END PARALLEL
  !
  dbldelta(1:nb,1:nb,1:nk_local) = dbldelta(1:nb,1:nb,1:nk_local) / DBLE(6 * nkBZ)
  !
END SUBROUTINE libtetrabz_dbldelta_main
!
! 2nd step of tetrahedra method.
!
SUBROUTINE libtetrabz_dbldelta2(nb,ej,w)
  !
  USE libtetrabz_common, ONLY : libtetrabz_sort
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: nb
  REAL(8),INTENT(IN) :: ej(3,nb)
  REAL(8),INTENT(INOUT) :: w(nb,3)
  !
  INTEGER :: ib, ii, indx(3)
  REAL(8) :: a(3,3), e(3), V
  character(len=500) :: msg
  !
  DO ib = 1, nb
     !
     IF(maxval(ABS(ej(1:3,ib))) < 1d-10) then
     ! MG reduce tolerance wrt original version.
     !IF(maxval(ABS(ej(1:3,ib))) < 1d-22) then
      !write(std_out, ej(1:3,ib))
      write(msg, *)"Nesting for band index:", ib, "ej:", ej(1:3,ib)
      MSG_WARNING(msg)
      !MSG_ERROR(msg)
      !MSG_ERROR("STOP Nesting !!")
     end if
     !
     w(ib,1:3) = 0d0
     e(1:3) = ej(1:3, ib)
     CALL libtetrabz_sort(3,e,indx)
     !
     DO ii = 1, 3
        a(1:3,ii) = (0d0 - e(ii)) / (e(1:3) - e(ii))
     END DO
     !
     IF((e(1) < 0d0 .AND. 0d0 <= e(2)) .OR. (e(1) <= 0d0 .AND. 0d0 < e(2))) THEN
        !
        !V = a(2,1) * a(3,1) / (0d0 - e(1))
        V = a(2,1)           / (e(3) - e(1))
        !
        w(ib,indx(1)) = V * (a(1,2) + a(1,3))
        w(ib,indx(2)) = V * a(2,1)
        w(ib,indx(3)) = V * a(3,1)
        !
     ELSE IF((e(2) <= 0d0 .AND. 0d0 < e(3)) .OR. (e(2) < 0d0 .AND. 0d0 <= e(3))) THEN
        !
        !V = a(1,3) * a(2,3) / (e(3) - 0d0)
        V = a(2,3)           / (e(3) - e(1))
        !
        w(ib,indx(1)) = V * a(1,3)
        w(ib,indx(2)) = V * a(2,3)
        w(ib,indx(3)) = V * (a(3,1) + a(3,2))
        !
     END IF
     !
  END DO ! ib
  !
END SUBROUTINE libtetrabz_dbldelta2
!
END MODULE libtetrabz_dbldelta_mod
