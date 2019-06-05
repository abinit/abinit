
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
MODULE libtetrabz_polstat_mod

  use defs_basis
  use m_errors
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: libtetrabz_polstat
  !
CONTAINS
!
! Compute Static polalization function
!
SUBROUTINE libtetrabz_polstat(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,comm) 
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
     CALL libtetrabz_polstat_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig1,eig2,nk_local,wghtd)
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
     CALL libtetrabz_polstat_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig1,eig2,nk_local,wght)
  END IF
  !
  ABI_FREE(ik_global)
  ABI_FREE(ik_local)
  !
END SUBROUTINE libtetrabz_polstat
!
! Main SUBROUTINE for polalization function : Theta(- E1) * Theta(E2) / (E2 - E1)
!
SUBROUTINE libtetrabz_polstat_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig1,eig2,nk_local,polstat)
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
  REAL(8),INTENT(IN) :: wlsm(4,20), eig1(nb,nkBZ), eig2(nb,nkBZ)
  REAL(8),INTENT(OUT) :: polstat(nb,nb,nk_local)
  !
  INTEGER :: ib, it, indx(4)
  REAL(8) :: e(4), ei1(4,nb), ei2(4), ej1(4,nb), ej2(4,nb), thr = 1d-10, tsmall(4,4), V, w1(nb,4), w2(nb,4)
  !
  polstat(1:nb,1:nb,1:nk_local) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(eig1,eig2,ik_global,ik_local,nb,nt_local,polstat,thr,wlsm) &
  !$OMP & PRIVATE(e,ei1,ei2,ej1,ej2,ib,indx,it,tsmall,V,w1,w2)
  !
  DO it = 1, nt_local
     !
     DO ib = 1, nb
        ei1(1:4, ib) = MATMUL(wlsm(1:4,1:20), eig1(ib, ik_global(1:20,it)))
        ej1(1:4, ib) = MATMUL(wlsm(1:4,1:20), eig2(ib, ik_global(1:20,it)))
     END DO
     !
     !$OMP DO
     DO ib = 1, nb
        !
        w1(1:nb,1:4) = 0d0
        e(1:4) = ei1(1:4, ib)
        CALL libtetrabz_sort(4,e,indx)
        !
        IF(e(1) <= 0d0 .AND. 0d0 < e(2)) THEN
           !
           CALL libtetrabz_tsmall_a1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polstat2(nb,ei2,ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF(e(2) <= 0d0 .AND. 0d0 < e(3)) THEN
           !
           CALL libtetrabz_tsmall_b1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polstat2(nb,ei2,ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b2(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polstat2(nb,ei2,ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b3(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polstat2(nb,ei2,ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF(e(3) <= 0d0 .AND. 0d0 < e(4)) THEN
           !
           CALL libtetrabz_tsmall_c1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polstat2(nb,ei2,ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_c2(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polstat2(nb,ei2,ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_c3(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polstat2(nb,ei2,ej2,w2)
              w1(1:nb,indx(1:4)) = w1(1:nb,            indx(1:4)) &
              &       + V * MATMUL(w2(1:nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF(e(4) <= 0d0) THEN
           !
           ei2(1:4     ) = ei1(1:4,  ib)
           ej2(1:4,1:nb) = ej1(1:4,1:nb)
           CALL libtetrabz_polstat2(nb,ei2,ej2,w2)
           w1(1:nb,1:4) = w1(1:nb,1:4) + w2(1:nb,1:4)
           !
        ELSE
           !
           CYCLE
           !
        END IF
        !
        polstat(1:nb,ib,ik_local(1:20,it)) = polstat(1:nb,ib,   ik_local(1:20,it)) &
        &                                + MATMUL(w1(1:nb,1:4), wlsm(1:4,1:20))
        !
     END DO ! ib
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !
  !$OMP END PARALLEL
  !
  polstat(1:nb,1:nb,1:nk_local) = polstat(1:nb,1:nb,1:nk_local) / DBLE(6 * nkBZ)
  !
END SUBROUTINE libtetrabz_polstat_main
!
! Tetrahedra method for theta( - E2)
!
SUBROUTINE libtetrabz_polstat2(nb,ei1,ej1,w1)
  !
  USE libtetrabz_common, ONLY : libtetrabz_sort, &
  &                             libtetrabz_tsmall_a1, libtetrabz_tsmall_b1, &
  &                             libtetrabz_tsmall_b2, libtetrabz_tsmall_b3, &
  &                             libtetrabz_tsmall_c1, libtetrabz_tsmall_c2, &
  &                             libtetrabz_tsmall_c3
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: nb
  REAL(8),INTENT(IN) :: ei1(4), ej1(4,nb)
  REAL(8),INTENT(INOUT) :: w1(nb,4)
  !
  INTEGER :: ib, indx(4)
  REAL(8) :: de(4), e(4), thr = 1d-8, tsmall(4,4), V, w2(4)
  !
  DO ib = 1, nb
     !
     w1(ib,1:4) = 0d0
     e(1:4) = - ej1(1:4, ib)
     CALL libtetrabz_sort(4,e,indx)
     !
     IF((e(1) <= 0d0 .AND. 0d0 < e(2)) .OR. (e(1) < 0d0 .AND. 0d0 <= e(2))) THEN
        !
        CALL libtetrabz_tsmall_a1(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib) - ei1(indx(1:4)))
           CALL libtetrabz_polstat3(de,w2)
           w1(ib,indx(1:4)) = w1(ib,                    indx(1:4)) &
           &                + V * MATMUL(w2(1:4), tsmall(1:4,1:4))
           !
        END IF
        !
     ELSE IF((e(2) <= 0d0 .AND. 0d0 < e(3)) .OR. (e(2) < 0d0 .AND. 0d0 <= e(3))) THEN
        !
        CALL libtetrabz_tsmall_b1(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib) - ei1(indx(1:4)))
           CALL libtetrabz_polstat3(de,w2)
           w1(ib,indx(1:4)) = w1(ib,indx(1:4)) &
           &          + V * MATMUL(w2(        1:4 ), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_b2(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib) - ei1(indx(1:4)))
           CALL libtetrabz_polstat3(de,w2)
           w1(ib,indx(1:4)) = w1(ib,                    indx(1:4)) &
           &                + V * MATMUL(w2(1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_b3(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib) - ei1(indx(1:4)))
           CALL libtetrabz_polstat3(de,w2)
           w1(ib,indx(1:4)) = w1(ib,                    indx(1:4)) &
           &                + V * MATMUL(w2(1:4), tsmall(1:4,1:4))
           !
        END IF
        !
     ELSE IF((e(3) <= 0d0 .AND. 0d0 < e(4)) .OR. (e(3) < 0d0 .AND. 0d0 <= e(4))) THEN
        !
        CALL libtetrabz_tsmall_c1(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib) - ei1(indx(1:4)))
           CALL libtetrabz_polstat3(de,w2)
           w1(ib,indx(1:4)) = w1(ib,                    indx(1:4)) &
           &                + V * MATMUL(w2(1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_c2(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib) - ei1(indx(1:4)))
           CALL libtetrabz_polstat3(de,w2)
           w1(ib,indx(1:4)) = w1(ib,                    indx(1:4)) &
           &                + V * MATMUL(w2(1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_c3(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib) - ei1(indx(1:4)))
           CALL libtetrabz_polstat3(de,w2)
           w1(ib,indx(1:4)) = w1(ib,                    indx(1:4)) &
           &                + V * MATMUL(w2(1:4), tsmall(1:4,1:4))
           !
        END IF
        !
     ELSE IF( e(4) <= 0d0 ) THEN
        !
        de(1:4) = ej1(1:4,ib) - ei1(1:4)
        CALL libtetrabz_polstat3(de,w2)
        w1(ib,1:4) = w1(ib,1:4) + w2(1:4)
        !
     END IF
     !
  END DO ! ib = 1, nb
  !
END SUBROUTINE libtetrabz_polstat2
!
! Tetarahedra method for delta(om - ep + e)
!
SUBROUTINE libtetrabz_polstat3(de,w1)
  !
  USE libtetrabz_common, ONLY : libtetrabz_sort
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: de(4)
  REAL(8),INTENT(INOUT) :: w1(4)
  !
  INTEGER :: ii, indx(4)
  REAL(8) :: e(4), ln(4), thr, thr2
  !
  e(1:4) = de(1:4)
  CALL libtetrabz_sort(4,e,indx)
  !
  thr = MAXVAL(e(1:4)) * 1d-3
  thr2 = 1d-8
  !
  DO ii = 1, 4
     IF(e(ii) < thr2) THEN
        IF(ii == 3) THEN
           STOP "  Nesting ! "
        END IF
        ln(ii) = 0d0
        e(ii) = 0d0
     ELSE
        ln(ii) = LOG(e(ii))
     END IF
  END DO
  !
  IF(ABS(e(4) - e(3)) < thr ) THEN
     IF(ABS(e(4) - e(2)) < thr ) THEN
        IF(ABS(e(4) - e(1)) < thr ) THEN
           !
           ! e(4) = e(3) = e(2) = e(1)
           !
           w1(indx(4)) = 0.25d0 / e(4)
           w1(indx(3)) = w1(indx(4))
           w1(indx(2)) = w1(indx(4))
           w1(indx(1)) = w1(indx(4))
           !
        ELSE
           !
           ! e(4) = e(3) = e(2)
           !
           w1(indx(4)) = libtetrabz_polstat_1211(e(4),e(1),ln(4),ln(1))
           w1(indx(3)) = w1(indx(4))
           w1(indx(2)) = w1(indx(4))
           w1(indx(1)) = libtetrabz_polstat_1222(e(1),e(4),ln(1),ln(4))
           !
           IF(ANY(w1(1:4) < 0d0)) THEN
              WRITE(std_out,'(100e15.5)') e(1:4)
              WRITE(std_out,'(100e15.5)') w1(indx(1:4))
              MSG_ERROR("weighting 4=3=2")
           END IF
           !
        END IF
     ELSE IF(ABS(e(2) - e(1)) < thr) THEN
        !
        ! e(4) = e(3), e(2) = e(1)
        !
        w1(indx(4)) = libtetrabz_polstat_1221(e(4),e(2), ln(4),ln(2))
        w1(indx(3)) = w1(indx(4))
        w1(indx(2)) = libtetrabz_polstat_1221(e(2),e(4), ln(2),ln(4))
        w1(indx(1)) = w1(indx(2))
        !
        IF(ANY(w1(1:4) < 0d0)) THEN
           WRITE(std_out,'(100e15.5)') e(1:4)
           WRITE(std_out,'(100e15.5)') w1(indx(1:4))
           MSG_ERROR("weighting 4=3 2=1")
        END IF
        !
     ELSE
        !
        ! e(4) = e(3)
        !
        w1(indx(4)) = libtetrabz_polstat_1231(e(4),e(1),e(2),ln(4),ln(1),ln(2))
        w1(indx(3)) = w1(indx(4))
        w1(indx(2)) = libtetrabz_polstat_1233(e(2),e(1),e(4),ln(2),ln(1),ln(4))
        w1(indx(1)) = libtetrabz_polstat_1233(e(1),e(2),e(4),ln(1),ln(2),ln(4))
        !
        IF(ANY(w1(1:4) < 0d0)) THEN
           WRITE(std_out,'(100e15.5)') e(1:4)
           WRITE(std_out,'(100e15.5)') w1(indx(1:4))
           MSG_ERROR("weighting 4=3")
        END IF
        !
     END IF
  ELSE IF(ABS(e(3) - e(2)) < thr) THEN
     IF(ABS(e(3) - e(1)) < thr) THEN
        !
        ! e(3) = e(2) = e(1)
        !
        w1(indx(4)) = libtetrabz_polstat_1222(e(4),e(3), ln(4),ln(3))
        w1(indx(3)) = libtetrabz_polstat_1211(e(3),e(4), ln(3),ln(4))
        w1(indx(2)) = w1(indx(3))
        w1(indx(1)) = w1(indx(3))
        !
        IF(ANY(w1(1:4) < 0d0)) THEN
           WRITE(std_out,'(100e15.5)') e(1:4)
           WRITE(std_out,'(100e15.5)') w1(indx(1:4))
           MSG_ERROR("weighting 3=2=1")
        END IF
        !
     ELSE
        !
        ! e(3) = e(2)
        !
        w1(indx(4)) = libtetrabz_polstat_1233(e(4),e(1),e(3),ln(4),ln(1),ln(3))
        w1(indx(3)) = libtetrabz_polstat_1231(e(3),e(1),e(4),ln(3),ln(1),ln(4))
        w1(indx(2)) = w1(indx(3))
        w1(indx(1)) = libtetrabz_polstat_1233(e(1),e(4),e(3),ln(1),ln(4),ln(3))
        !
        IF(ANY(w1(1:4) < 0d0)) THEN
           WRITE(std_out,'(100e15.5)') e(1:4)
           WRITE(std_out,'(100e15.5)') w1(indx(1:4))
           MSG_ERROR("weighting 3=2")
        END IF
        !
     END IF
  ELSE IF(ABS(e(2) - e(1)) < thr) THEN
     !
     ! e(2) = e(1)
     !
     w1(indx(4)) = libtetrabz_polstat_1233(e(4),e(3),e(2),ln(4),ln(3),ln(2))
     w1(indx(3)) = libtetrabz_polstat_1233(e(3),e(4),e(2),ln(3),ln(4),ln(2))
     w1(indx(2)) = libtetrabz_polstat_1231(e(2),e(3),e(4),ln(2),ln(3),ln(4))
     w1(indx(1)) = w1(indx(2))
     !
     IF(ANY(w1(1:4) < 0d0)) THEN
        WRITE(std_out,'(100e15.5)') e(1:4)
        WRITE(std_out,'(100e15.5)') w1(indx(1:4))
        MSG_ERROR("weighting 2=1")
     END IF
     !
  ELSE
     !
     ! Different each other.
     !
     w1(indx(4)) = libtetrabz_polstat_1234(e(4),e(1),e(2),e(3),ln(4),ln(1),ln(2),ln(3))
     w1(indx(3)) = libtetrabz_polstat_1234(e(3),e(1),e(2),e(4),ln(3),ln(1),ln(2),ln(4))
     w1(indx(2)) = libtetrabz_polstat_1234(e(2),e(1),e(3),e(4),ln(2),ln(1),ln(3),ln(4))
     w1(indx(1)) = libtetrabz_polstat_1234(e(1),e(2),e(3),e(4),ln(1),ln(2),ln(3),ln(4))
     !
     IF(ANY(w1(1:4) < 0d0)) THEN
        WRITE(std_out,'(100e15.5)') e(1:4)
        WRITE(std_out,'(100e15.5)') w1(indx(1:4))
        MSG_ERROR("weighting")
     END IF
     !
  END IF
  !
END SUBROUTINE libtetrabz_polstat3
!
! Results of Integration (1-x-y-z)/(g0+(g1-g0)x+(g2-g0)y+(g3-g0))
!  for 0<x<1, 0<y<1-x, 0<z<1-x-y
!
! 1, Different each other
!
FUNCTION libtetrabz_polstat_1234(g1,g2,g3,g4,lng1,lng2,lng3,lng4)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1,g2,g3,g4,lng1,lng2,lng3,lng4
  REAL(8) :: w, libtetrabz_polstat_1234
  !
  REAL(8) :: w2, w3, w4
  !
  w2 = ((lng2 - lng1)/(g2 - g1)*g2 - 1d0)*g2/(g2 - g1)
  w3 = ((lng3 - lng1)/(g3 - g1)*g3 - 1d0)*g3/(g3 - g1)
  w4 = ((lng4 - lng1)/(g4 - g1)*g4 - 1d0)*g4/(g4 - g1)
  w2 = ((w2 - w3)*g2)/(g2 - g3)
  w4 = ((w4 - w3)*g4)/(g4 - g3)
  w = (w4 - w2)/(g4 - g2)
  libtetrabz_polstat_1234 = w
  !
END FUNCTION libtetrabz_polstat_1234
!
! 2, g4 = g1
!
FUNCTION libtetrabz_polstat_1231(g1,g2,g3,lng1,lng2,lng3)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1,g2,g3,lng1,lng2,lng3
  REAL(8) :: w, libtetrabz_polstat_1231
  !
  REAL(8) :: w2, w3
  !
  w2 = ((lng2 - lng1)/(g2 - g1)*g2 - 1d0)*g2**2/(g2 - g1) - g1/( &
  &   2d0)
  w2 = w2/(g2 - g1)
  w3 = ((lng3 - lng1)/(g3 - g1)*g3 - 1d0)*g3**2/(g3 - g1) - g1/( &
  &   2d0)
  w3 = w3/(g3 - g1)
  w = (w3 - w2)/(g3 - g2)
  libtetrabz_polstat_1231 = w
  !
END FUNCTION libtetrabz_polstat_1231
!
! 3, g4 = g3
!
FUNCTION libtetrabz_polstat_1233(g1,g2,g3,lng1,lng2,lng3)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1,g2,g3,lng1,lng2,lng3
  REAL(8) :: w,libtetrabz_polstat_1233
  !
  REAL(8) :: w2, w3
  !
  w2 = (lng2 - lng1)/(g2 - g1)*g2 - 1d0
  w2 = (g2*w2)/(g2 - g1)
  w3 = (lng3 - lng1)/(g3 - g1)*g3 - 1d0
  w3 = (g3*w3)/(g3 - g1)
  w2 = (w3 - w2)/(g3 - g2)
  w3 = (lng3 - lng1)/(g3 - g1)*g3 - 1d0
  w3 = 1d0 - (2d0*w3*g1)/(g3 - g1)
  w3 = w3/(g3 - g1)
  w = (g3*w3 - g2*w2)/(g3 - g2)
  libtetrabz_polstat_1233 = w
  !
END FUNCTION libtetrabz_polstat_1233
!
! 4, g4 = g1 and g3 = g2
!
FUNCTION libtetrabz_polstat_1221(g1,g2,lng1,lng2)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1, g2, lng1, lng2
  REAL(8) :: w, libtetrabz_polstat_1221
  !
  w = 1d0 - (lng2 - lng1)/(g2 - g1)*g1
  w = -1d0 + (2d0*g2*w)/(g2 - g1)
  w = -1d0 + (3d0*g2*w)/(g2 - g1)
  w = w/(2d0*(g2 - g1))
  libtetrabz_polstat_1221 = w
  !
END FUNCTION libtetrabz_polstat_1221
!
! 5, g4 = g3 = g2
!
FUNCTION libtetrabz_polstat_1222(g1,g2,lng1,lng2)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1, g2, lng1, lng2
  REAL(8) :: w, libtetrabz_polstat_1222
  !
  w = (lng2 - lng1)/(g2 - g1)*g2 - 1d0
  w = (2d0*g1*w)/(g2 - g1) - 1d0
  w = (3d0*g1*w)/(g2 - g1) + 1d0
  w = w/(2d0*(g2 - g1))
  libtetrabz_polstat_1222 = w
  !
END FUNCTION libtetrabz_polstat_1222
!
! 6, g4 = g3 = g1
!
FUNCTION libtetrabz_polstat_1211(g1,g2,lng1,lng2)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1,g2,lng1,lng2
  REAL(8) :: w, libtetrabz_polstat_1211
  !
  w = -1d0 + (lng2 - lng1)/(g2 - g1)*g2
  w = -1d0 + (2d0*g2*w)/(g2 - g1)
  w = -1d0 + (3d0*g2*w)/(2d0*(g2 - g1))
  w = w/(3d0*(g2 - g1))
  libtetrabz_polstat_1211 = w
  !
END FUNCTION libtetrabz_polstat_1211
!
END MODULE libtetrabz_polstat_mod
