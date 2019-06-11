
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
MODULE libtetrabz_fermigr_mod
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: libtetrabz_fermigr
  !
CONTAINS
!
! Compute Fermi's golden rule
!
SUBROUTINE libtetrabz_fermigr(ltetra,bvec,nb,nge,eig1,eig2,ngw,wght,ne,e0,comm) 
  !
  USE ISO_C_BINDING
  USE libtetrabz_common, ONLY : libtetrabz_initialize, libtetrabz_interpol_indx, libtetrabz_mpisum_dv
  IMPLICIT NONE
  !
  INTEGER(C_INT),INTENT(IN) :: ltetra, nb, nge(3), ngw(3), ne
  REAL(C_DOUBLE),INTENT(IN) :: bvec(9), e0(ne), eig1(nb,PRODUCT(nge(1:3))), eig2(nb,PRODUCT(nge(1:3)))
  REAL(C_DOUBLE),INTENT(OUT) :: wght(ne*nb*nb,PRODUCT(ngw(1:3)))
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
     ABI_MALLOC(wghtd, (ne*nb*nb,1,nk_local))
     CALL libtetrabz_fermigr_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig1,eig2,ne,e0,nk_local,wghtd)
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
     IF(PRESENT(comm)) CALL libtetrabz_mpisum_dv(comm, ne*nb*nb*PRODUCT(ngw(1:3)), wght)
     !
  ELSE
     CALL libtetrabz_fermigr_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig1,eig2,ne,e0,nk_local,wght)
  END IF
  !
  ABI_FREE(ik_global)
  ABI_FREE(ik_local)
  !
END SUBROUTINE libtetrabz_fermigr
!
! Main SUBROUTINE for Fermi's Gorlden rule : Theta(- E1) * Theta(E2) * Delta(E2 - E1 - w)
!
SUBROUTINE libtetrabz_fermigr_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig1,eig2,ne,e0,nk_local,fermigr)
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
  REAL(8),INTENT(IN) :: wlsm(4,20), eig1(nb,nkBZ), eig2(nb,nkBZ), e0(ne)
  REAL(8),INTENT(OUT) :: fermigr(ne*nb,nb,nk_local)
  !
  INTEGER :: it, ib, indx(4)
  REAL(8) :: e(4), ei1(4,nb), ei2(4), ej1(4,nb), ej2(4,nb), thr = 1d-10, tsmall(4,4), V, w1(ne*nb,4), w2(ne*nb,4)
  !
  fermigr(1:ne*nb,1:nb,1:nk_local) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(eig1,eig2,e0,fermigr,ik_global,ik_local,nb,ne,nt_local,thr,wlsm) &
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
              CALL libtetrabz_fermigr2(nb,ne,e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,     1:4 ), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF( e(2) <= 0d0 .AND. 0d0 < e(3)) THEN
           !
           CALL libtetrabz_tsmall_b1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_fermigr2(nb,ne,e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,     1:4 ), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b2(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_fermigr2(nb,ne,e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,     1:4 ), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_b3(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_fermigr2(nb,ne,e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,     1:4 ), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF( e(3) <= 0d0 .AND. 0d0 < e(4)) THEN
           !
           CALL libtetrabz_tsmall_c1(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_fermigr2(nb,ne,e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,     1:4 ), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_c2(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_fermigr2(nb,ne,e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,     1:4 ), tsmall(1:4,1:4))
              !
           END IF
           !
           CALL libtetrabz_tsmall_c3(e,V,tsmall)
           !
           IF(V > thr) THEN
              !
              ei2(1:4     ) = MATMUL(tsmall(1:4,1:4), ei1(indx(1:4),  ib))
              ej2(1:4,1:nb) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),1:nb))
              CALL libtetrabz_fermigr2(nb,ne,e0,ei2,ej2,w2)
              w1(1:ne*nb,indx(1:4)) = w1(1:ne*nb,indx(1:4)) &
              &          + V * MATMUL(w2(1:ne*nb,     1:4 ), tsmall(1:4,1:4))
              !
           END IF
           !
        ELSE IF(e(4) <= 0d0) THEN
           !
           ei2(1:4     ) = ei1(1:4,  ib)
           ej2(1:4,1:nb) = ej1(1:4,1:nb)
           CALL libtetrabz_fermigr2(nb,ne,e0,ei2,ej2,w2)
           w1(1:ne*nb,1:4) = w1(1:ne*nb,1:4) + w2(1:ne*nb,1:4)
           !
        END IF
        !
        fermigr(1:ne*nb,ib,ik_local(1:20,it)) = fermigr(1:ne*nb,ib,   ik_local(1:20,it)) &
        &                                   + MATMUL(w1(1:ne*nb,1:4), wlsm(1:4,1:20))
        !
     END DO ! ib = 1, nb
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !
  !$OMP END PARALLEL
  !
  fermigr(1:ne*nb,1:nb,1:nk_local) = fermigr(1:ne*nb,1:nb,1:nk_local) / DBLE(6 * nkBZ)
  !
END SUBROUTINE libtetrabz_fermigr_main
!
! Tetrahedra method for theta( - E2)
!
SUBROUTINE libtetrabz_fermigr2(nb,ne,e0,ei1,ej1,w1)
  !
  USE libtetrabz_common, ONLY : libtetrabz_sort, &
  &                             libtetrabz_tsmall_a1, libtetrabz_tsmall_b1, &
  &                             libtetrabz_tsmall_b2, libtetrabz_tsmall_b3, &
  &                             libtetrabz_tsmall_c1, libtetrabz_tsmall_c2, &
  &                             libtetrabz_tsmall_c3
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: nb, ne
  REAL(8),INTENT(IN) :: e0(ne), ei1(4), ej1(4,nb)
  REAL(8),INTENT(OUT) :: w1(ne,nb,4)
  !
  INTEGER :: ib, indx(4)
  REAL(8) :: de(4), e(4), thr = 1d-8, tsmall(4,4), V, w2(ne,4)
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
           CALL libtetrabz_fermigr3(ne,e0,de,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,        1:4 ), tsmall(1:4,1:4))
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
           CALL libtetrabz_fermigr3(ne,e0,de,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,        1:4 ), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_b2(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib) - ei1(indx(1:4)))
           CALL libtetrabz_fermigr3(ne,e0,de,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,        1:4 ), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_b3(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib) - ei1(indx(1:4)))
           CALL libtetrabz_fermigr3(ne,e0,de,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,        1:4 ), tsmall(1:4,1:4))
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
           CALL libtetrabz_fermigr3(ne,e0,de,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,        1:4 ), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_c2(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib) - ei1(indx(1:4)))
           CALL libtetrabz_fermigr3(ne,e0,de,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,        1:4 ), tsmall(1:4,1:4))
           !
        END IF
        !
        CALL libtetrabz_tsmall_c3(e,V,tsmall)
        !
        IF(V > thr) THEN
           !
           de(1:4) = MATMUL(tsmall(1:4,1:4), ej1(indx(1:4),ib) - ei1(indx(1:4)))
           CALL libtetrabz_fermigr3(ne,e0,de,w2)
           w1(1:ne,ib,indx(1:4)) = w1(1:ne,ib,indx(1:4)) &
           &          + V * MATMUL(w2(1:ne,        1:4 ), tsmall(1:4,1:4))
           !
        END IF
        !
     ELSE IF(e(4) <= 0d0) THEN
        !
        de(1:4) = ej1(1:4,ib) - ei1(1:4)
        CALL libtetrabz_fermigr3(ne,e0,de,w2)
        w1(1:ne,ib,1:4) = w1(1:ne,ib,1:4) + w2(1:ne,1:4)
        !
     END IF
     !
  END DO ! ib = 1, nb
  !
END SUBROUTINE libtetrabz_fermigr2
!
!
!
SUBROUTINE libtetrabz_fermigr3(ne,e0,de,w1)
  !
  USE libtetrabz_common, ONLY : libtetrabz_sort, &
  &                             libtetrabz_triangle_a1, libtetrabz_triangle_b1, &
  &                             libtetrabz_triangle_b2, libtetrabz_triangle_c1
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: ne
  REAL(8),INTENT(IN) :: e0(ne), de(4)
  REAL(8),INTENT(OUT) :: w1(ne,4)
  !
  INTEGER :: ie, indx(4)
  REAL(8) :: e(4), tsmall(3,4), V, w2(3)
  !
  w2(1:3) = 1d0 / 3d0
  !
  w1(1:ne,1:4) = 0d0
  e(1:4) = de(1:4)
  CALL libtetrabz_sort(4,e,indx)
  !
  DO ie = 1, ne
     !
     IF(e(1) < e0(ie) .AND. e0(ie) <= e(2)) THEN
        !
        CALL libtetrabz_triangle_a1(e(1:4) - e0(ie),V,tsmall)
        w1(ie,indx(1:4)) = w1(ie,indx(1:4)) + V * SUM(tsmall(1:3,1:4), 1) / 3d0
        !
     ELSE IF(e(2) < e0(ie) .AND. e0(ie) <= e(3)) THEN
        !
        CALL libtetrabz_triangle_b1(e(1:4) - e0(ie),V,tsmall)
        w1(ie,indx(1:4)) = w1(ie,indx(1:4)) + V * SUM(tsmall(1:3,1:4), 1) / 3d0
        !
        CALL libtetrabz_triangle_b2(e(1:4) - e0(ie),V,tsmall)
        w1(ie,indx(1:4)) = w1(ie,indx(1:4)) + V * SUM(tsmall(1:3,1:4), 1) / 3d0
        !
     ELSE IF(e(3) < e0(ie) .AND. e0(ie) < e(4)) THEN
        !
        CALL libtetrabz_triangle_c1(e(1:4) - e0(ie),V,tsmall)
        w1(ie,indx(1:4)) = w1(ie,indx(1:4)) + V * SUM(tsmall(1:3,1:4), 1) / 3d0
        !
     END IF
     !
  END DO ! ie
  !
END SUBROUTINE libtetrabz_fermigr3
!
END MODULE libtetrabz_fermigr_mod
