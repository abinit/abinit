#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_FFTHyboffdiag
!! NAME
!!  m_FFTHyboffdiag
!! 
!! FUNCTION 
!!  Almost useless. Just uses for FFT time evolution
!!  of number of electrons
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

#define FFTHyboffdiag_FORWARD   1
#define FFTHyboffdiag_BACKWARD -1
#include "defs.h"
MODULE m_FFTHyboffdiag
USE m_global

IMPLICIT NONE

!!***

!!****t* m_FFTHyboffdiag/FFTHyboffdiag
!! NAME
!!  FFTHyboffdiag
!!
!! FUNCTION
!!  This structured datatype contains the necessary data
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

TYPE FFTHyboffdiag
  LOGICAL :: set = .FALSE.
  INTEGER :: size
  DOUBLE PRECISION :: Ts
  DOUBLE PRECISION :: fs
  INTEGER         , ALLOCATABLE, DIMENSION(:) :: bit_rev   
  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:) :: data_inout
END TYPE FFTHyboffdiag
!!***

CONTAINS
!!***

!!****f* ABINIT/m_FFTHyboffdiag/FFTHyboffdiag_init
!! NAME
!!  FFTHyboffdiag_init
!!
!! FUNCTION
!!  Initialize ...
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=FFT
!!  n=Number of point (should be power of 2)
!!  samples_sec=number of samples per sec
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE FFTHyboffdiag_init(op,n,samples_sec)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'FFTHyboffdiag_init'
!End of the abilint section

  TYPE(FFTHyboffdiag), INTENT(INOUT) :: op
  INTEGER     , INTENT(IN   ) :: n
  DOUBLE PRECISION, INTENT(IN   ) :: samples_sec
!Local variables ------------------------------
  INTEGER :: i
  INTEGER :: inv_bit
  INTEGER :: total_size

  ! Check n is a power of 2
  IF ( n .LT. 2 .OR. IAND(n,n-1) .NE. 0 ) THEN
    CALL WARNALL("FFTHyboffdiag_init : array size is not a power of 2 -> auto fix")
    i = 1
    DO WHILE ( i .LT. n )
      i = ISHFT( i, 1 )
    END DO
    total_size = ISHFT(i, -1)
  ELSE
    total_size = n
  END IF
  
  op%size = total_size
  op%Ts = DBLE(total_size) / samples_sec
  op%fs = 1.d0 / DBLE(total_size)

  FREEIF(op%data_inout)
  MALLOC(op%data_inout,(0:total_size-1))
  op%data_inout(0:total_size-1) = CMPLX(0.d0,0.d0,8)
  FREEIF(op%bit_rev)
  MALLOC(op%bit_rev,(0:total_size-1))
  op%bit_rev = 0

  DO i = 1, total_size-1
    inv_bit = FFTHyboffdiag_mirror(i,total_size)
    op%bit_rev(inv_bit) = i 
  END DO

  op%set = .TRUE.

END SUBROUTINE FFTHyboffdiag_init
!!***

!!****f* ABINIT/m_FFTHyboffdiag/FFTHyboffdiag_mirror
!! NAME
!!  FFTHyboffdiag_mirror
!!
!! FUNCTION
!!  mirror bits of an integer
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  i=bits
!!  n=integer
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

INTEGER FUNCTION FFTHyboffdiag_mirror(i,n)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'FFTHyboffdiag_mirror'
!End of the abilint section

  INTEGER, INTENT(IN ) :: i
  INTEGER, INTENT(IN ) :: n
!Local variable -------------------------------
  INTEGER              :: icp
  INTEGER              :: ncp

  icp = i
  ncp = n

  FFTHyboffdiag_mirror = 0

  DO WHILE( ncp .GT. 1 )
    FFTHyboffdiag_mirror = IOR(ISHFT(FFTHyboffdiag_mirror,1) , IAND(icp,1))
    icp = ISHFT(icp,-1)
    ncp = ISHFT(ncp,-1)
  END DO

END FUNCTION FFTHyboffdiag_mirror
!!***

!!****f* ABINIT/m_FFTHyboffdiag/FFTHyboffdiag_setData
!! NAME
!!  FFTHyboffdiag_setData
!!
!! FUNCTION
!!  set input data (in time)
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=FFThyb
!!  array_in=data in
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE FFTHyboffdiag_setData(op, array_in)
!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'FFTHyboffdiag_setData'
!End of the abilint section

  TYPE(FFTHyboffdiag), INTENT(INOUT) :: op
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: array_in
  INTEGER :: size_in
  INTEGER :: i

  size_in = SIZE(array_in)
  op%data_inout = CMPLX(0.d0,0.d0)
!  IF ( size_in .NE. op%size ) &
!    CALL WARNALL("FFTHyboffdiag_setData : size_in != size")
  
  DO i = 0, MIN(op%size,size_in)-1
    op%data_inout(i) = CMPLX(array_in(i+1), 0.d0)
  END DO

END SUBROUTINE FFTHyboffdiag_setData
!!***

!!****f* ABINIT/m_FFTHyboffdiag/FFTHyboffdiag_run
!! NAME
!!  FFTHyboffdiag_run
!!
!! FUNCTION
!!  perform FFT
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=FFT
!!  dir=direction FFTHyboffdiag_FORWARD or FFTHyboffdiag_BACKWARD
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE FFTHyboffdiag_run(op, dir)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'FFTHyboffdiag_run'
!End of the abilint section

  TYPE(FFTHyboffdiag), INTENT(INOUT) :: op
  INTEGER     , INTENT(IN   ) :: dir
!Local variables ------------------------------
  INTEGER :: imax
  INTEGER :: istep
  INTEGER :: i
  INTEGER :: m
  INTEGER :: j
  DOUBLE PRECISION :: wtmp
  DOUBLE PRECISION :: wr
  DOUBLE PRECISION :: wpr
  DOUBLE PRECISION :: wpi
  DOUBLE PRECISION :: wi
  DOUBLE PRECISION :: theta
  DOUBLE PRECISION :: twoPi
  COMPLEX(KIND=8) :: tc
 
  imax = 1;
  istep = 2;
  
  twoPi = DBLE(dir)*2.d0*ACOS(-1.d0)
 
  DO WHILE ( imax .LT. op%size )
    istep = ISHFT(imax,1)
    theta = twoPi/DBLE(istep)
    Wtmp  = SIN(0.5d0*theta)
    wpr   = -2.d0*wtmp*wtmp
    wpi   = SIN(theta)
    wr    = 1.0d0
    wi    = 0.d0
    DO m = 0, imax-1
      DO i = m, op%size-1, istep
        j= i+imax
        tc = CMPLX( wr* REAL (op%data_inout(op%bit_rev(j)))  &
                   -wi* AIMAG(op%data_inout(op%bit_rev(j))), &
                    wr* AIMAG(op%data_inout(op%bit_rev(j)))  &
                   +wi* REAL (op%data_inout(op%bit_rev(j))), 8 )
        op%data_inout(op%bit_rev(j)) = op%data_inout(op%bit_rev(i)) - tc
        op%data_inout(op%bit_rev(i)) = op%data_inout(op%bit_rev(i)) + tc
      END DO
      wtmp = wr
      wr   = wr*wpr - wi*wpi  +wr
      wi   = wi*wpr + wtmp*wpi+wi
    END DO
    imax = istep
  END DO
  IF ( dir .EQ. FFTHyboffdiag_FORWARD ) &
    op%data_inout = op%data_inout*op%fs

END SUBROUTINE FFTHyboffdiag_run
!!***

!!****f* ABINIT/m_FFTHyboffdiag/FFTHyboffdiag_getData
!! NAME
!!  FFTHyboffdiag_getData
!!
!! FUNCTION
!!  get result
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=FFT
!!
!! OUTPUT
!!  bound=number of frequencies
!!  array_out=output data
!!  freqs=output frequencies
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE FFTHyboffdiag_getData(op, bound, array_out, freqs)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'FFTHyboffdiag_getData'
!End of the abilint section

  TYPE(FFTHyboffdiag), INTENT(IN) :: op
  INTEGER         ,               INTENT(OUT) :: bound
  DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: array_out
  DOUBLE PRECISION, DIMENSION(:), OPTIONAL, INTENT(OUT) :: freqs
!Local variables ------------------------------
  INTEGER :: i
  INTEGER :: size_out
  DOUBLE PRECISION :: re
  DOUBLE PRECISION :: im

  size_out = SIZE(array_out)
!  IF ( SIZE(array_out) .NE. op%size/2 ) &
!    CALL WARNALL("FFTHyboffdiag_getData : size_in != size")
  array_out = 0.d0
  bound = MIN(op%size/2, size_out)
  DO i=1, bound
    re = REAL(op%data_inout(op%bit_rev(i-1)))
    im = AIMAG(op%data_inout(op%bit_rev(i-1)))
    array_out(i) = SQRT(re*re + im*im)
  END DO
  IF ( PRESENT( freqs ) .AND. bound .LE. SIZE(freqs)) THEN
    DO i=1, bound
        freqs(i) = DBLE(i-1)/op%Ts 
    END DO
!  ELSE IF ( PRESENT( freqs ) .AND. bound .GT. SIZE(freqs) ) THEN
!    CALL WARNALL("FFHyb_getData : freqs does is too small")
  END IF
  
END SUBROUTINE FFTHyboffdiag_getData
!!***

!!****f* ABINIT/m_FFTHyboffdiag/FFTHyboffdiag_destroy
!! NAME
!!  FFTHyboffdiag_destroy
!!
!! FUNCTION
!!  destroy every thing
!!
!! COPYRIGHT
!!  Copyright (C) 2013 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  op=FFT
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE FFTHyboffdiag_destroy(op)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'FFTHyboffdiag_destroy'
!End of the abilint section

  TYPE(FFTHyboffdiag), INTENT(INOUT) :: op

  FREEIF(op%bit_rev)
  FREEIF(op%data_inout)
  op%set = .FALSE.
END SUBROUTINE FFTHyboffdiag_destroy
!!***

END MODULE m_FFTHyboffdiag
!!***

