
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
MODULE libtetrabz_polcmplx_mod
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: libtetrabz_polcmplx
  !
CONTAINS
!
! Compute Polarization of imaginary frequency
!
SUBROUTINE libtetrabz_polcmplx(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,ne,e0,comm) 
  !
  USE ISO_C_BINDING
  USE libtetrabz_common, ONLY : libtetrabz_initialize, libtetrabz_interpol_indx, libtetrabz_mpisum_zv
  IMPLICIT NONE
  !
  INTEGER(C_INT),INTENT(IN) :: ltetra, nb, nge(3), ngw(3), ne
  REAL(C_DOUBLE),INTENT(IN) :: bvec(9), eig1(nb,PRODUCT(nge(1:3))), eig2(nb,PRODUCT(nge(1:3)))
  COMPLEX(C_DOUBLE_COMPLEX),INTENT(IN) :: e0(ne)
  COMPLEX(C_DOUBLE_COMPLEX),INTENT(OUT) :: wght(ne*nb*nb,PRODUCT(ngw(1:3)))
  INTEGER(C_INT),INTENT(IN),OPTIONAL :: comm
  !
  LOGICAL :: linterpol
  INTEGER :: nt_local, nk_local, nkBZ, ik, kintp(20), nintp
  INTEGER,ALLOCATABLE :: ik_global(:,:), ik_local(:,:)
  REAL(8) :: wlsm(4,20), wintp(1,20)
  REAL(8),ALLOCATABLE :: kvec(:,:)
  COMPLEX(8),ALLOCATABLE :: wghtd(:,:,:)
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
     ABI_MALLOC(wghtd, (ne*nb*nb,1,nk_local))
     CALL libtetrabz_polcmplx_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig1,eig2,ne,e0,nk_local,wghtd)
     !
     ! Interpolation
     !
     wght(1:ne*nb*nb,1:PRODUCT(ngw(1:3))) = 0d0
     DO ik = 1, nk_local
        CALL libtetrabz_interpol_indx(nintp,ngw,kvec(1:3,ik),kintp,wintp)
        wght(1:ne*nb*nb,kintp(1:nintp)) = wght(1:ne*nb*nb,             kintp(1:nintp)) &
        &                       + MATMUL(wghtd(1:ne*nb*nb,1:1,ik), wintp(1:1,1:nintp))
     END DO ! ik = 1, nk_local
     ABI_FREE(wghtd)
     ABI_FREE(kvec)
     !
     IF(PRESENT(comm)) CALL libtetrabz_mpisum_zv(comm, ne*nb*nb*PRODUCT(ngw(1:3)), wght)
     !
  ELSE
     CALL libtetrabz_polcmplx_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig1,eig2,ne,e0,nk_local,wght)
  END IF
  !
  ABI_FREE(ik_global)
  ABI_FREE(ik_local)
  !
END SUBROUTINE libtetrabz_polcmplx
!
! Main SUBROUTINE for Polaization (Imaginaly axis) : Theta(- E1) * Theta(E2) / (E2 - E1 - iw)
!
SUBROUTINE libtetrabz_polcmplx_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig1,eig2,ne,e0,nk_local,polcmplx)
  !
  USE libtetrabz_common, ONLY : libtetrabz_sort, &
  &                             libtetrabz_tsmall_a1, libtetrabz_tsmall_b1, &
  &                             libtetrabz_tsmall_b2, libtetrabz_tsmall_b3, &
  &                             libtetrabz_tsmall_c1, libtetrabz_tsmall_c2, &
  &                             libtetrabz_tsmall_c3
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: nt_local, nb, nkBZ, nk_local, ne, &
  &                     ik_global(20,nt_local), ik_local(20,nt_local)
  REAL(8),INTENT(IN) :: wlsm(4,20), eig1(nb,nkBZ), eig2(nb,nkBZ)
  COMPLEX(8),INTENT(IN) :: e0(ne)
  COMPLEX(8),INTENT(OUT) :: polcmplx(ne*nb,nb,nk_local)
  !
  INTEGER :: ib, indx(4), it
  REAL(8) :: e(4), ei1(4,nb), ej1(4,nb), ei2(4), ej2(4,nb), thr = 1d-8, tsmall(4,4), V
  COMPLEX(8) :: w1(ne*nb,4), w2(ne*nb,4)
  !
  polcmplx(1:ne*nb,1:nb,1:nk_local) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(eig1,eig2,e0,ik_global,ik_local,nb,ne,nt_local,polcmplx,thr,wlsm) &
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
        w1(1:ne*nb,1:4) = 0d0
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
              CALL libtetrabz_polcmplx2(nb,ne,e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,            indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
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
              CALL libtetrabz_polcmplx2(nb,ne,e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,            indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b2(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polcmplx2(nb,ne,e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,            indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b3(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polcmplx2(nb,ne,e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,            indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
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
              CALL libtetrabz_polcmplx2(nb,ne,e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,            indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_c2(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polcmplx2(nb,ne,e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,            indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_c3(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_polcmplx2(nb,ne,e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,            indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,1:4), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF( e(4) <= 0d0 ) THEN
           !
           ei2(1:4     ) = ei1(1:4,  ib)
           ej2(1:4,1:nb) = ej1(1:4,1:nb)
           CALL libtetrabz_polcmplx2(nb,ne,e0,ei2,ej2,w2)
           w1(1:ne*nb,1:4) = w1(1:ne*nb,1:4) + w2(1:ne*nb,1:4)
           !
        ELSE
           !
           CYCLE
           !
        END IF
        !
        polcmplx(1:ne*nb,ib,ik_local(1:20,it)) = polcmplx(1:ne*nb,ib,   ik_local(1:20,it)) &
        &                                     + MATMUL(w1(1:ne*nb,1:4), wlsm(1:4,1:20))
        !
     END DO ! ib = 1, nb
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !
  !$OMP END PARALLEL
  !
  polcmplx(1:ne*nb,1:nb,1:nk_local) = polcmplx(1:ne*nb,1:nb,1:nk_local) / DBLE(6 * nkBZ)
  !
END SUBROUTINE libtetrabz_polcmplx_main
!
! Tetrahedra method for theta( - E2)
!
SUBROUTINE libtetrabz_polcmplx2(nb,ne,e0,ei1,ej1,w1)
  !
  USE libtetrabz_common, ONLY : libtetrabz_sort, &
  &                             libtetrabz_tsmall_a1, libtetrabz_tsmall_b1, &
  &                             libtetrabz_tsmall_b2, libtetrabz_tsmall_b3, &
  &                             libtetrabz_tsmall_c1, libtetrabz_tsmall_c2, &
  &                             libtetrabz_tsmall_c3
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: nb, ne
  COMPLEX(8),INTENT(IN) :: e0(ne)
  REAL(8),INTENT(IN) :: ei1(4), ej1(4,nb)
  COMPLEX(8),INTENT(OUT) :: w1(ne,nb,4)
  !
  INTEGER :: ib, indx(4)
  REAL(8) :: de(4), e(4), thr = 1d-8, tsmall(4,4), V
  COMPLEX(8) :: w2(ne,4)
  !
  DO ib = 1, nb
     !
     w1(1:ne,ib,1:4) = 0d0
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
           CALL libtetrabz_polcmplx3(ne,e0,de,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,         indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
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
           CALL libtetrabz_polcmplx3(ne,e0,de,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,         indx(1:4)) &
           &            + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_b2(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib) - ei1(indx(1:4)))
           CALL libtetrabz_polcmplx3(ne,e0,de,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,         indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_b3(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib) - ei1(indx(1:4)))
           CALL libtetrabz_polcmplx3(ne,e0,de,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,         indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
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
           CALL libtetrabz_polcmplx3(ne,e0,de,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,         indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_c2(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib) - ei1(indx(1:4)))
           CALL libtetrabz_polcmplx3(ne,e0,de,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,         indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_c3(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib) - ei1(indx(1:4)))
           CALL libtetrabz_polcmplx3(ne,e0,de,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,         indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,1:4), tsmall(1:4,1:4))
           !
        END IF
        !
     ELSE IF(e(4) <= 0d0) THEN
        !
        de(1:4) = ej1(1:4,ib) - ei1(1:4)
        CALL libtetrabz_polcmplx3(ne,e0,de,w2)
        w1(1:ne,ib,1:4) = w1(1:ne,ib,1:4) + w2(1:ne,1:4)
        !
     END IF
     !
  END DO
  !
END SUBROUTINE libtetrabz_polcmplx2
!
! Tetarahedra method for delta(om - ep + e)
!
SUBROUTINE libtetrabz_polcmplx3(ne,e0,de,w1)
  !
  USE libtetrabz_common, ONLY : libtetrabz_sort
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: ne
  COMPLEX(8),INTENT(IN) :: e0(ne)
  REAL(8),INTENT(IN) :: de(4)
  COMPLEX(8),INTENT(OUT) :: w1(ne,4)
  !
  INTEGER :: ie, indx(4)
  REAL(8) :: e(4), thr, w2(2,4), x(4)
  !
  e(1:4) = de(1:4)
  CALL libtetrabz_sort(4,e,indx)
  !
  DO ie = 1, ne
     !
     ! I don't know which one is better.
     ! The former is more stable.
     ! The latter is more accurate ?
     !
     w1(ie,1:4) = 0.25d0 / (de(1:4) + e0(ie))
     !
     CYCLE
     !
     x(1:4) = (e(1:4) + DBLE(e0(ie))) / AIMAG(e0(ie))
     !thr = maxval(de(1:4)) * 1d-3
     thr = max(1d-3,  MAXVAL(x(1:4)) * 1d-2)
     !
     IF(ABS(x(4) - x(3)) < thr) THEN
        IF(ABS(x(4) - x(2)) < thr) THEN
           IF(ABS(x(4) - x(1)) < thr) THEN
              !
              ! e(4) = e(3) = e(2) = e(1)
              !
              w2(1,4) = 0.25d0 * x(4) / ((1d0 + x(4)**2))
              w2(2,4) = 0.25d0        / ((1d0 + x(4)**2))
              w2(1:2,3) = w2(1:2,4)
              w2(1:2,2) = w2(1:2,4)
              w2(1:2,1) = w2(1:2,4)
              !
           ELSE
              !
              ! e(4) = e(3) = e(2)
              !
              w2(1:2,4) = libtetrabz_polcmplx_1211(x(4),x(1))
              w2(1:2,3) = w2(1:2,4)
              w2(1:2,2) = w2(1:2,4)
              w2(1:2,1) = libtetrabz_polcmplx_1222(x(1),x(4))
              !
              !IF(ANY(w2(1:2,1:4) < 0d0)) THEN
              !   WRITE(*,*) ie
              !   WRITE(*,'(100e15.5)') x(1:4)
              !   WRITE(*,'(2e15.5)') w2(1:2,1:4)
              !   STOP "weighting 4=3=2"
              !END IF
              !
           END IF
        ELSE IF(ABS(x(2) - x(1)) < thr ) THEN
           !
           ! e(4) = e(3), e(2) = e(1)
           !
           w2(1:2,4) = libtetrabz_polcmplx_1221(x(4),x(2))
           w2(1:2,3) = w2(1:2,4)
           w2(1:2,2) = libtetrabz_polcmplx_1221(x(2),x(4))
           w2(1:2,1) = w2(1:2,2)
           !
           !IF(ANY(w2(1:2,1:4) < 0d0)) THEN
           !   WRITE(*,*) ie
           !   WRITE(*,'(100e15.5)') x(1:4)
           !   WRITE(*,'(2e15.5)') w2(1:2,1:4)
           !   STOP "weighting 4=3 2=1"
           !END IF
           !
        ELSE
           !
           ! e(4) = e(3)
           !
           w2(1:2,4) = libtetrabz_polcmplx_1231(x(4),x(1),x(2))
           w2(1:2,3) = w2(1:2,4)
           w2(1:2,2) = libtetrabz_polcmplx_1233(x(2),x(1),x(4))
           w2(1:2,1) = libtetrabz_polcmplx_1233(x(1),x(2),x(4))
           !
           !IF(ANY(w2(1:2,1:4) < 0d0)) THEN
           !   WRITE(*,*) ie
           !   WRITE(*,'(100e15.5)') x(1:4)
           !   WRITE(*,'(2e15.5)') w2(1:2,1:4)
           !   STOP "weighting 4=3"
           !END IF
           !
        END IF
     ELSE IF(ABS(x(3) - x(2)) < thr) THEN
        IF(ABS(x(3) - x(1)) < thr) THEN
           !
           ! e(3) = e(2) = e(1)
           !
           w2(1:2,4) = libtetrabz_polcmplx_1222(x(4),x(3))
           w2(1:2,3) = libtetrabz_polcmplx_1211(x(3),x(4))
           w2(1:2,2) = w2(1:2,3)
           w2(1:2,1) = w2(1:2,3)
           !
           !IF(ANY(w2(1:2,1:4) < 0d0)) THEN
           !   WRITE(*,*) ie
           !   WRITE(*,'(100e15.5)') x(1:4)
           !   WRITE(*,'(2e15.5)') w2(1:2,1:4)
           !   STOP "weighting 3=2=1"
           !END IF
           !
        ELSE
           !
           ! e(3) = e(2)
           !
           w2(1:2,4) = libtetrabz_polcmplx_1233(x(4),x(1),x(3))
           w2(1:2,3) = libtetrabz_polcmplx_1231(x(3),x(1),x(4))
           w2(1:2,2) = w2(1:2,3)
           w2(1:2,1) = libtetrabz_polcmplx_1233(x(1),x(4),x(3))
           !
           !IF(ANY(w2(1:2,1:4) < 0d0)) THEN
           !   WRITE(*,*) ie
           !   WRITE(*,'(100e15.5)') x(1:4)
           !   WRITE(*,'(2e15.5)') w2(1:2,1:4)
           !   STOP "weighting 3=2"
           !END IF
           !
        END IF
     ELSE IF(ABS(x(2) - x(1)) < thr) THEN
        !
        ! e(2) = e(1)
        !
        w2(1:2,4) = libtetrabz_polcmplx_1233(x(4),x(3),x(2))
        w2(1:2,3) = libtetrabz_polcmplx_1233(x(3),x(4),x(2))
        w2(1:2,2) = libtetrabz_polcmplx_1231(x(2),x(3),x(4))
        w2(1:2,1) = w2(1:2,2)
        !
        !IF(ANY(w2(1:2,1:4) < 0d0)) THEN
        !   WRITE(*,*) ie
        !   WRITE(*,'(100e15.5)') x(1:4)
        !   WRITE(*,'(2e15.5)') w2(1:2,1:4)
        !   STOP "weighting 2=1"
        !END IF
        !
     ELSE
        !
        ! Different each other.
        !
        w2(1:2,4) = libtetrabz_polcmplx_1234(x(4),x(1),x(2),x(3))
        w2(1:2,3) = libtetrabz_polcmplx_1234(x(3),x(1),x(2),x(4))
        w2(1:2,2) = libtetrabz_polcmplx_1234(x(2),x(1),x(3),x(4))
        w2(1:2,1) = libtetrabz_polcmplx_1234(x(1),x(2),x(3),x(4))
        !
        !IF(ANY(w2(1:2,1:4) < 0d0)) THEN
        !   WRITE(*,*) ie
        !   WRITE(*,'(100e15.5)') x(1:4)
        !   WRITE(*,'(2e15.5)') w2(1:2,1:4)
        !   STOP "weighting"
        !END IF
        !
     END IF
     !
     w1(ie,indx(1:4)) = CMPLX(w2(1,1:4) /    AIMAG(e0(ie)), &
     &                        w2(2,1:4) / (- AIMAG(e0(ie))), KIND(0d0))
     !
  END DO ! ie
  !
END SUBROUTINE libtetrabz_polcmplx3
!
! Results of Integration (1-x-y-z)/(g0+(g1-g0)x+(g2-g0)y+(g3-g0))
!  for 0<x<1, 0<y<1-x, 0<z<1-x-y
!
! 1, Different each other
!
FUNCTION libtetrabz_polcmplx_1234(g1,g2,g3,g4)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1, g2, g3, g4
  REAL(8) :: w(2), libtetrabz_polcmplx_1234(2)
  !
  REAL(8) :: w2, w3, w4
  !
  ! Real
  !
  w2 = 2d0*(3d0*g2**2 - 1d0)*(ATAN(g2) - ATAN(g1)) + (g2**2 - &
  &      3d0)*g2*LOG((1d0 + g2**2)/( 1d0 + g1**2))
  w2 = -2d0*(g2**2 - 1d0) + w2/(g2 - g1 )
  w2 = w2/(g2 - g1 )
  w3 = 2d0*(3d0*g3**2 - 1d0)*(ATAN(g3) - ATAN(g1)) + (g3**2 -  &
  &      3d0)*g3*LOG((1d0 + g3**2)/( 1d0 + g1**2))
  w3 = -2d0*(g3**2 - 1d0) + w3/(g3 - g1 )
  w3 = w3/(g3 - g1 )
  w4 = 2d0*(3d0*g4**2 - 1d0)*(ATAN(g4) - ATAN(g1)) + (g4**2 -  &
  &      3d0)*g4*LOG((1d0 + g4**2)/( 1d0 + g1**2))
  w4 = -2d0*(g4**2 - 1d0) + w4/(g4 - g1 )
  w4 = w4/(g4 - g1 )
  w2 = (w2 - w3)/(g2 - g3)
  w4 = (w4 - w3)/(g4 - g3)
  w(1) = (w4 - w2)/(2d0*(g4 - g2))
  !
  ! Imaginal
  !
  w2 = 2d0*(3d0 - g2**2)* &
  &    g2*(ATAN(g2) - ATAN(g1)) + (3d0*g2**2 - 1d0)* &
  &    LOG((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 4d0*g2 - w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0*(3d0 - g3**2)* &
  &    g3*(ATAN(g3) - ATAN(g1)) + (3d0*g3**2 - 1d0)* &
  &    LOG((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 4d0*g3 - w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w4 = 2d0*(3d0 - g4**2)* &
  &    g4*(ATAN(g4) - ATAN(g1)) + (3d0*g4**2 - 1d0)* &
  &    LOG((1d0 + g4**2)/(1d0 + g1**2))
  w4 = 4d0*g4 - w4/(g4 - g1)
  w4 = w4/(g4 - g1)
  w2 = (w2 - w3)/(g2 - g3)
  w4 = (w4 - w3)/(g4 - g3)
  w(2) = (w4 - w2)/(2d0*(g4 - g2))
  !
  libtetrabz_polcmplx_1234 = w
END FUNCTION libtetrabz_polcmplx_1234
!
! 2, g4 = g1
!
FUNCTION libtetrabz_polcmplx_1231(g1,g2,g3) 
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1, g2, g3
  REAL(8) :: w(2), libtetrabz_polcmplx_1231(2)
  !
  REAL(8) :: w2, w3
  !
  ! Real
  !
  w2 = 2d0*(-1d0 + 3d0*g2**2)*(ATAN(g2) - ATAN(g1)) +  &
  &   g2*(-3d0 + g2**2)*LOG((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 2d0*(1d0 - g2**2) + w2/(g2 - g1)
  w2 = -g1 + w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0*(-1d0 + 3d0*g3**2)*(ATAN(g3) - ATAN(g1)) +  &
  &   g3*(-3d0 + g3**2)*LOG((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 2d0*(1 - g3**2) + w3/(g3 - g1)
  w3 = -g1 + w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w(1) = (w3 - w2)/(2d0*(g3 - g2))
  !
  ! Imaginal
  !
  w2 = 2d0* &
  &    g2*(3d0 - g2**2)*(ATAN(g2) - ATAN(g1)) + (-1d0 + 3d0*g2**2)* &
  &    LOG((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 4d0*g2 - w2/(g2 - g1)
  w2 = 1 + w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0* &
  &    g3*(3d0 - g3**2)*(ATAN(g3) - ATAN(g1)) + (-1d0 + 3d0*g3**2)* &
  &    LOG((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 4d0*g3 - w3/(g3 - g1)
  w3 = 1 + w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w(2) = (w3 - w2)/(2d0*(g3 - g2))

  libtetrabz_polcmplx_1231 = w
  !
END FUNCTION libtetrabz_polcmplx_1231
!
! 3, g4 = g3
!
FUNCTION libtetrabz_polcmplx_1233(g1, g2, g3)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1, g2, g3
  REAL(8) :: w(2), libtetrabz_polcmplx_1233(2)
  !
  REAL(8) :: w2, w3
  !
  ! Real
  !
  w2 = 2d0*(1d0 - 3d0*g2**2)*(ATAN(g2) - ATAN(g1)) +  &
  &   g2*(3d0 - g2**2)*LOG((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 2d0*(1 - g2**2) - w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0*(1d0 - 3d0*g3**2)*(ATAN(g3) - ATAN(g1)) +  &
  &   g3*(3d0 - g3**2)*LOG((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 2d0*(1 - g3**2) - w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w2 = (w3 - w2)/(g3 - g2)
  w3 = 4d0*(1d0 - 3d0*g1*g3)*(ATAN(g3) - ATAN(g1)) + (3d0*g1 +  &
  &      3d0*g3 - 3d0*g1*g3**2 + g3**3) * LOG((1d0 + g3**2)/( &
  &     1d0 + g1**2))
  w3 = -4d0*(1d0 - g1**2) + w3/(g3 - g1)
  w3 = 4d0*g1 + w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w(1) = (w3 - w2)/(2d0*(g3 - g2))
  !
  ! Imaginal
  !
  w2 = 2d0* &
  &    g2*(3d0 - g2**2)*(ATAN(g2) - ATAN(g1)) + (-1d0 + 3d0*g2**2)* &
  &    LOG((1d0 + g2**2)/(1d0 + g1**2))
  w2 = 4d0*g2 - w2/(g2 - g1)
  w2 = w2/(g2 - g1)
  w3 = 2d0* &
  &    g3*(3d0 - g3**2)*(ATAN(g3) - ATAN(g1)) + (-1d0 + 3d0*g3**2)* &
  &    LOG((1d0 + g3**2)/(1d0 + g1**2))
  w3 = 4d0*g3 - w3/(g3 - g1)
  w3 = w3/(g3 - g1)
  w2 = (w3 - w2)/(g3 - g2)
  w3 = (3d0*g1 - 3d0*g1*g3**2 + 3d0*g3 + g3**3)*(ATAN(g3) -  &
  &      ATAN(g1)) + (3d0*g1*g3 - 1d0)* &
  &    LOG((1d0 + g3**2)/(1d0 + g1**2))
  w3 = w3/(g3 - g1) - 4d0*g1
  w3 = w3/(g3 - g1) - 2d0
  w3 = (2d0*w3)/(g3 - g1)
  w(2) = (w3 - w2)/(2d0*(g3 - g2))
  !
  libtetrabz_polcmplx_1233 = w
END FUNCTION libtetrabz_polcmplx_1233
!
! 4, g4 = g1 and g3 = g2
!
FUNCTION libtetrabz_polcmplx_1221(g1,g2)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1, g2
  REAL(8) :: w(2),libtetrabz_polcmplx_1221(2)
  !
  ! Real
  !
  w(1) = -2d0*(-1d0 + 2d0*g1*g2 + g2**2)*(ATAN(g2) -  &
  &      ATAN(g1)) + (g1 + 2d0*g2 - g1*g2**2)* &
  &    LOG((1d0 + g2**2)/(1d0 + g1**2))
  w(1) = 2d0*(-1d0 + g1**2) + w(1)/(g2 - g1)
  w(1) = 3d0*g1 + w(1)/(g2 - g1)
  w(1) = 2d0 + (3d0*w(1))/(g2 - g1)
  w(1) = w(1)/(2d0*(g2 - g1))
  !
  ! Imaginal
  !
  w(2) = 2d0*(g1 + 2d0*g2 - g1*g2**2)*(ATAN(g2) -  &
  &      ATAN(g1)) + (-1d0 + 2d0*g1*g2 + g2**2)* &
  &    LOG((1 + g2**2)/(1 + g1**2))
  w(2) = -4d0*g1 + w(2)/(g2 - g1)
  w(2) = -3d0 + w(2)/(g2 - g1)
  w(2) = (3d0*w(2))/(2d0*(g2 - g1)**2)

  libtetrabz_polcmplx_1221 = w
  !
END FUNCTION libtetrabz_polcmplx_1221
!
! 5, g4 = g3 = g2
!
FUNCTION libtetrabz_polcmplx_1222(g1,g2)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1, g2
  REAL(8) :: w(2), libtetrabz_polcmplx_1222(2)
  !
  ! Real
  !
  w(1) = 2d0*(-1d0 + g1**2 + 2d0*g1*g2)*(ATAN(g2) -  &
  &      ATAN(g1)) + (-2d0*g1 - g2 + g1**2*g2) * LOG((1d0 + g2**2)/( &
  &     1d0 + g1**2))
  w(1) = 2d0*(1d0 - g1**2) + w(1)/(g2 - g1)
  w(1) = g1 - w(1)/(g2 - g1)
  w(1) = 1d0 - (3d0*w(1))/(g2 - g1)
  w(1) = w(1)/(2d0*(g2 - g1))
  !
  ! Imaginal
  !
  w(2) = 2d0*(-2d0*g1 - g2 + g1**2*g2)*(ATAN(g2) - ATAN(g1)) + (1d0 - &
  &       g1**2 - 2d0*g1*g2) * LOG((1d0 + g2**2)/(1d0 + g1**2))
  w(2) = 4d0*g1 + w(2)/(g2 - g1)
  w(2) = 1d0 + w(2)/(g2 - g1)
  w(2) = (3d0*w(2))/(2d0*(g2 - g1)**2)
  !
  libtetrabz_polcmplx_1222 = w
END FUNCTION libtetrabz_polcmplx_1222
!
! 6, g4 = g3 = g1
!
FUNCTION libtetrabz_polcmplx_1211(g1,g2)
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: g1, g2
  REAL(8) :: w(2), libtetrabz_polcmplx_1211(2)
  !
  ! Real
  !
  w(1) = 2d0*(3d0*g2**2 - 1d0)*(ATAN(g2) - ATAN(g1)) +  &
  &   g2*(g2**2 - 3d0)*LOG((1d0 + g2**2)/(1d0 + g1**2))
  w(1) = 2d0*(1d0 - g1**2) + w(1)/(g2 - g1)
  w(1) = -5d0*g1 + w(1)/(g2 - g1)
  w(1) = -11d0 + (3d0*w(1))/(g2 - g1)
  w(1) = w(1)/(6d0*(g2 - g1))
  !
  ! Imaginal
  !
  w(2) = 2d0*g2*(-3d0 + g2**2)*(ATAN(g2) - ATAN(g1)) + (1d0 -  &
  &      3d0*g2**2)*LOG((1d0 + g2**2)/(1d0 + g1**2))
  w(2) = 4d0*g2 + w(2)/(g2 - g1)
  w(2) = 1d0 + w(2)/(g2 - g1)
  w(2) = w(2)/(2d0*(g2 - g1)**2)
  libtetrabz_polcmplx_1211 = w
  !
END FUNCTION libtetrabz_polcmplx_1211
!
END MODULE libtetrabz_polcmplx_mod
