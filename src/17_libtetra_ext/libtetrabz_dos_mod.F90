
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
MODULE libtetrabz_dos_mod
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: libtetrabz_dos, libtetrabz_intdos
  !
CONTAINS
!
! Compute DOS
!
SUBROUTINE libtetrabz_dos(ltetra,bvec,nb,nge,eig,ngw,wght,ne,e0,comm)
  !
  USE ISO_C_BINDING
  USE libtetrabz_common, ONLY : libtetrabz_initialize, libtetrabz_interpol_indx, libtetrabz_mpisum_dv
  IMPLICIT NONE
  !
  INTEGER(C_INT),INTENT(IN) :: ltetra, nb, nge(3), ngw(3), ne
  REAL(C_DOUBLE),INTENT(IN) :: bvec(9), eig(nb,PRODUCT(nge(1:3))), e0(ne)
  REAL(C_DOUBLE),INTENT(OUT) :: wght(ne*nb,PRODUCT(ngw(1:3)))
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
     ABI_MALLOC(wghtd, (ne*nb,1,nk_local))
     CALL libtetrabz_dos_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig,ne,e0,nk_local,wghtd)
     !
     ! Interpolation
     !
     wght(1:ne*nb,1:PRODUCT(ngw(1:3))) = 0d0
     DO ik = 1, nk_local
        CALL libtetrabz_interpol_indx(nintp,ngw,kvec(1:3,ik),kintp,wintp)
        wght(1:ne*nb,kintp(1:nintp)) = wght(1:ne*nb,             kintp(1:nintp)) &
        &                    + MATMUL(wghtd(1:ne*nb,1:1,ik), wintp(1:1,1:nintp))
     END DO ! ik = 1, nk_local
     ABI_FREE(wghtd)
     ABI_FREE(kvec)
     !
     IF(PRESENT(comm)) CALL libtetrabz_mpisum_dv(comm, ne*nb*PRODUCT(ngw(1:3)), wght)
     !
  ELSE
     CALL libtetrabz_dos_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig,ne,e0,nk_local,wght)
  END IF
  !
  ABI_FREE(ik_global)
  ABI_FREE(ik_local)
  !
END SUBROUTINE libtetrabz_dos
!
! Compute Integrated DOS
!
SUBROUTINE libtetrabz_intdos(ltetra,bvec,nb,nge,eig,ngw,wght,ne,e0,comm)
  !
  USE ISO_C_BINDING
  USE libtetrabz_common, ONLY : libtetrabz_initialize, libtetrabz_interpol_indx, libtetrabz_mpisum_dv
  IMPLICIT NONE
  !
  INTEGER(C_INT),INTENT(IN) :: ltetra, nb, nge(3), ngw(3), ne
  REAL(C_DOUBLE),INTENT(IN) :: bvec(9), eig(nb,PRODUCT(nge(1:3))), e0(ne)
  REAL(C_DOUBLE),INTENT(OUT) :: wght(ne*nb,PRODUCT(ngw(1:3)))
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
     ABI_MALLOC(wghtd, (ne*nb,1,nk_local))
     CALL libtetrabz_intdos_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig,ne,e0,nk_local,wghtd)
     !
     ! Interpolation
     !
     wght(1:ne*nb,1:PRODUCT(ngw(1:3))) = 0d0
     DO ik = 1, nk_local
        CALL libtetrabz_interpol_indx(nintp,ngw,kvec(1:3,ik),kintp,wintp)
        wght(1:ne*nb,kintp(1:nintp)) = wght(1:ne*nb,         kintp(1:nintp)) &
        &                    + MATMUL(wghtd(1:ne*nb,1:1,ik), wintp(1:1,1:nintp))
     END DO ! ik = 1, nk_local
     ABI_FREE(wghtd)
     ABI_FREE(kvec)
     !
     IF(PRESENT(comm)) CALL libtetrabz_mpisum_dv(comm, ne*nb*PRODUCT(ngw(1:3)), wght)
     !
  ELSE
     CALL libtetrabz_intdos_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig,ne,e0,nk_local,wght)
  END IF
  !
  ABI_FREE(ik_global)
  ABI_FREE(ik_local)
  !
END SUBROUTINE libtetrabz_intdos
!
! Main SUBROUTINE for Dos : Delta(E - E1)
!
SUBROUTINE libtetrabz_dos_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig,ne,e0,nk_local,dos)
  !
  USE libtetrabz_common, ONLY : libtetrabz_sort, &
  &                             libtetrabz_triangle_a1, libtetrabz_triangle_b1, &
  &                             libtetrabz_triangle_b2, libtetrabz_triangle_c1
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: nt_local, nb, nkBZ, nk_local, ne, &
  &                     ik_global(20,nt_local), ik_local(20,nt_local)
  REAL(8),INTENT(IN) :: wlsm(4,20), eig(nb,nkBZ), e0(ne)
  REAL(8),INTENT(OUT) :: dos(ne,nb,nk_local)
  !
  INTEGER :: ib, it, ie, indx(4)  !, ii
  REAL(8) :: e(4), ei1(4,nb), tsmall(3,4), V, w1(ne,4) !, wmul(ne,20)
  !
  dos(1:ne, 1:nb, 1:nk_local) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(dos,eig,e0,ik_global,ik_local,nb,ne,nt_local,wlsm) &
  !$OMP & PRIVATE(e,ei1,ib,ie,indx,it,tsmall,V,w1)
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
        w1(1:ne,1:4) = 0d0
        e(1:4) = ei1(1:4, ib)
        CALL libtetrabz_sort(4,e,indx)
        !
        DO ie = 1, ne
           !
           IF((e(1) <  e0(ie) .AND. e0(ie) <= e(2)) .OR. &
           &  (e(1) <= e0(ie) .AND. e0(ie) <  e(2))) THEN
              !
              CALL libtetrabz_triangle_a1(e(1:4) - e0(ie),V,tsmall)
              w1(ie,indx(1:4)) = w1(ie,indx(1:4)) + V * SUM(tsmall(1:3,1:4), 1) / 3d0
              !
           ELSE IF((e(2) <  e0(ie) .AND. e0(ie) <= e(3)) .OR. &
           &       (e(2) <= e0(ie) .AND. e0(ie) <  e(3))) THEN
              !
              CALL libtetrabz_triangle_b1(e(1:4) - e0(ie),V,tsmall)
              w1(ie,indx(1:4)) = w1(ie,indx(1:4)) + V * SUM(tsmall(1:3,1:4), 1) / 3d0
              !
              CALL libtetrabz_triangle_b2(e(1:4) - e0(ie),V,tsmall)
              w1(ie,indx(1:4)) = w1(ie,indx(1:4)) + V * SUM(tsmall(1:3,1:4), 1) / 3d0
              !
           ELSE IF((e(3) < e0(ie) .AND. e0(ie) <= e(4)) .OR. &
           &       (e(3) <= e0(ie) .AND. e0(ie) < e(4))) THEN
              !
              CALL libtetrabz_triangle_c1(e(1:4) - e0(ie),V,tsmall)
              w1(ie,indx(1:4)) = w1(ie,indx(1:4)) + V * SUM(tsmall(1:3,1:4), 1) / 3d0
              !
           END IF
           !
        END DO ! ie
        !

        !call DGEMM("N", "N", ne, 20, 4, 1d0, w1, ne, wlsm, 4, 0d0, wmul, ne)
        !wmul(:, 1:4) = MATMUL(w1, wlsm(:, 1:4))
        !do ii=1,4
        !    ie = ik_local(ii,it)
        !    dos(:,ib, ie) = dos(:,ib, ie) + wmul(:, ii)
        !end do

        ! This call is expensive and 1:20 can be reduced to 1:4 if Blochl
        dos(1:ne,ib,ik_local(1:20,it)) = dos(1:ne,ib,    ik_local(1:20,it)) &
        &                        + MATMUL(w1(1:ne, 1:4), wlsm(1:4,1:20))
        !
     END DO ! ib
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !
  !$OMP END PARALLEL
  !
  dos(1:ne,1:nb,1:nk_local) = dos(1:ne,1:nb,1:nk_local) / DBLE(6 * nkBZ)
  !
END SUBROUTINE libtetrabz_dos_main
!
! Main SUBROUTINE for integrated Dos : theta(E - E1)
!
SUBROUTINE libtetrabz_intdos_main(wlsm,nt_local,ik_global,ik_local,nb,nkBZ,eig,ne,e0,nk_local,intdos)
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
  REAL(8),INTENT(IN) :: wlsm(4,20), eig(nb,nkBZ), e0(ne)
  REAL(8),INTENT(OUT) :: intdos(ne,nb,nk_local)
  !
  INTEGER :: ib, it, ie, indx(4)
  REAL(8) :: e(4), ei1(4,nb), tsmall(4,4), V, w1(ne,4)
  !
  intdos(1:ne, 1:nb, 1:nk_local) = 0d0
  !
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & SHARED(eig,e0,ik_global,ik_local,intdos,nb,ne,nt_local,wlsm) &
  !$OMP & PRIVATE(e,ei1,ib,ie,indx,it,tsmall,V,w1)
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
        w1(1:ne,1:4) = 0d0
        e(1:4) = ei1(1:4, ib)
        CALL libtetrabz_sort(4,e,indx)
        !
        DO ie = 1, ne
           !
           IF((e(1) <= e0(ie) .AND. e0(ie) <  e(2)) .OR. &
           &  (e(1) <  e0(ie) .AND. e0(ie) <= e(2))) THEN
              !
              CALL libtetrabz_tsmall_a1(e - e0(ie),V,tsmall)
              w1(ie,indx(1:4)) = w1(ie,indx(1:4)) + V * SUM(tsmall(1:4,1:4), 1) * 0.25d0
              !
           ELSE IF((e(2) <= e0(ie) .AND. e0(ie) <  e(3)) .OR. &
           &       (e(2) <  e0(ie) .AND. e0(ie) <= e(3))) THEN
              !
              CALL libtetrabz_tsmall_b1(e - e0(ie),V,tsmall)
              w1(ie,indx(1:4)) = w1(ie,indx(1:4)) + V * SUM(tsmall(1:4,1:4), 1) * 0.25d0
              !
              CALL libtetrabz_tsmall_b2(e - e0(ie),V,tsmall)
              w1(ie,indx(1:4)) = w1(ie,indx(1:4)) + V * SUM(tsmall(1:4,1:4), 1) * 0.25d0
              !
              CALL libtetrabz_tsmall_b3(e - e0(ie),V,tsmall)
              w1(ie,indx(1:4)) = w1(ie,indx(1:4)) + V * SUM(tsmall(1:4,1:4), 1) * 0.25d0
              !
           ELSE IF((e(3) <= e0(ie) .AND. e0(ie) <  e(4)) .OR. &
           &       (e(3) <  e0(ie) .AND. e0(ie) <= e(4))) THEN
              !
              CALL libtetrabz_tsmall_c1(e - e0(ie),V,tsmall)
              w1(ie,indx(1:4)) = w1(ie,indx(1:4)) + V * SUM(tsmall(1:4,1:4), 1) * 0.25d0
              !
              CALL libtetrabz_tsmall_c2(e - e0(ie),V,tsmall)
              w1(ie,indx(1:4)) = w1(ie,indx(1:4)) + V * SUM(tsmall(1:4,1:4), 1) * 0.25d0
              !
              CALL libtetrabz_tsmall_c3(e - e0(ie),V,tsmall)
              w1(ie,indx(1:4)) = w1(ie,indx(1:4)) + V * SUM(tsmall(1:4,1:4), 1) * 0.25d0
              !
           ELSE IF(e(4) <= e0(ie)) THEN
              !
              w1(ie,1:4) = 0.25d0
              !
           ELSE
              !
              CYCLE
              !
           END IF
           !
        END DO ! ie
        !
        intdos(1:ne,ib,ik_local(1:20,it)) = intdos(1:ne,ib,    ik_local(1:20,it)) &
        &                              + MATMUL(w1(1:ne, 1:4), wlsm(1:4,1:20))
        !
     END DO ! ib
     !$OMP END DO NOWAIT
     !
  END DO ! it
  !
  !$OMP END PARALLEL
  !
  intdos(1:ne,1:nb,1:nk_local) = intdos(1:ne,1:nb,1:nk_local) / DBLE(6 * nkBZ)
  !
END SUBROUTINE libtetrabz_intdos_main
!
END MODULE libtetrabz_dos_mod
