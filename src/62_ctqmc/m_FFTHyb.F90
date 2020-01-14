
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!!****m* ABINIT/m_FFTHyb
!! NAME
!!  m_FFTHyb
!! 
!! FUNCTION 
!!  Almost useless. Just uses for FFT time evolution
!!  of number of electrons
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
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

#define FFTHyb_FORWARD   1
#define FFTHyb_BACKWARD -1
#include "defs.h"
MODULE m_FFTHyb
USE m_global

IMPLICIT NONE

!!***

PRIVATE

!!****t* m_FFTHyb/FFTHyb
!! NAME
!!  FFTHyb
!!
!! FUNCTION
!!  This structured datatype contains the necessary data
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

TYPE, PUBLIC :: FFTHyb
  LOGICAL _PRIVATE :: set = .FALSE.
  INTEGER          :: size
  DOUBLE PRECISION _PRIVATE :: Ts
  DOUBLE PRECISION _PRIVATE :: fs
  INTEGER         , ALLOCATABLE, DIMENSION(:) _PRIVATE :: bit_rev   
  COMPLEX(KIND=8) , ALLOCATABLE, DIMENSION(:) _PRIVATE :: data_inout
END TYPE FFTHyb
!!***

PUBLIC  :: FFTHyb_init
PRIVATE :: FFTHyb_mirror
PUBLIC  :: FFTHyb_setData
PUBLIC  :: FFTHyb_run
PUBLIC  :: FFTHyb_getData
PUBLIC  :: FFTHyb_destroy


CONTAINS
!!***

!!****f* ABINIT/m_FFTHyb/FFTHyb_init
!! NAME
!!  FFTHyb_init
!!
!! FUNCTION
!!  Initialize ...
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=FFT
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

SUBROUTINE FFTHyb_init(this,n,samples_sec)

!Arguments ------------------------------------
  TYPE(FFTHyb), INTENT(INOUT) :: this
  INTEGER     , INTENT(IN   ) :: n
  DOUBLE PRECISION, INTENT(IN   ) :: samples_sec
!Local variables ------------------------------
  INTEGER :: i
  INTEGER :: inv_bit
  INTEGER :: total_size

  ! Check n is a power of 2
  IF ( n .LT. 2 .OR. IAND(n,n-1) .NE. 0 ) THEN
    CALL WARNALL("FFTHyb_init : array size is not a power of 2 -> auto fix")
    i = 1
    DO WHILE ( i .LT. n )
      i = ISHFT( i, 1 )
    END DO
    total_size = ISHFT(i, -1)
  ELSE
    total_size = n
  END IF
  
  this%size = total_size
  this%Ts = DBLE(total_size) / samples_sec
  this%fs = 1.d0 / DBLE(total_size)

  FREEIF(this%data_inout)
  MALLOC(this%data_inout,(0:total_size-1))
  this%data_inout(0:total_size-1) = CMPLX(0.d0,0.d0,8)
  FREEIF(this%bit_rev)
  MALLOC(this%bit_rev,(0:total_size-1))
  this%bit_rev = 0

  DO i = 1, total_size-1
    inv_bit = FFTHyb_mirror(i,total_size)
    this%bit_rev(inv_bit) = i 
  END DO

  this%set = .TRUE.

END SUBROUTINE FFTHyb_init
!!***

!!****f* ABINIT/m_FFTHyb/FFTHyb_mirror
!! NAME
!!  FFTHyb_mirror
!!
!! FUNCTION
!!  mirror bits of an integer
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
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

INTEGER FUNCTION FFTHyb_mirror(i,n)

!Arguments ------------------------------------
  INTEGER, INTENT(IN ) :: i
  INTEGER, INTENT(IN ) :: n
!Local variable -------------------------------
  INTEGER              :: icp
  INTEGER              :: ncp

  icp = i
  ncp = n

  FFTHyb_mirror = 0

  DO WHILE( ncp .GT. 1 )
    FFTHyb_mirror = IOR(ISHFT(FFTHyb_mirror,1) , IAND(icp,1))
    icp = ISHFT(icp,-1)
    ncp = ISHFT(ncp,-1)
  END DO

END FUNCTION FFTHyb_mirror
!!***

!!****f* ABINIT/m_FFTHyb/FFTHyb_setData
!! NAME
!!  FFTHyb_setData
!!
!! FUNCTION
!!  set input data (in time)
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=FFThyb
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

SUBROUTINE FFTHyb_setData(this, array_in)
!Arguments ------------------------------------
  TYPE(FFTHyb), INTENT(INOUT) :: this
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: array_in
  INTEGER :: size_in
  INTEGER :: i

  size_in = SIZE(array_in)
  this%data_inout = CMPLX(0.d0,0.d0,KIND=4)
!  IF ( size_in .NE. this%size ) &
!    CALL WARNALL("FFTHyb_setData : size_in != size")
  
  DO i = 0, MIN(this%size,size_in)-1
    this%data_inout(i) = CMPLX(array_in(i+1), 0.d0,8)
  END DO

END SUBROUTINE FFTHyb_setData
!!***

!!****f* ABINIT/m_FFTHyb/FFTHyb_run
!! NAME
!!  FFTHyb_run
!!
!! FUNCTION
!!  perform FFT
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=FFT
!!  dir=direction FFTHyb_FORWARD or FFTHyb_BACKWARD
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

SUBROUTINE FFTHyb_run(this, dir)

!Arguments ------------------------------------
  TYPE(FFTHyb), INTENT(INOUT) :: this
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
 
  DO WHILE ( imax .LT. this%size )
    istep = ISHFT(imax,1)
    theta = twoPi/DBLE(istep)
    Wtmp  = SIN(0.5d0*theta)
    wpr   = -2.d0*wtmp*wtmp
    wpi   = SIN(theta)
    wr    = 1.0d0
    wi    = 0.d0
    DO m = 0, imax-1
      DO i = m, this%size-1, istep
        j= i+imax
        tc = CMPLX( wr* REAL (this%data_inout(this%bit_rev(j)))  &
                   -wi* AIMAG(this%data_inout(this%bit_rev(j))), &
                    wr* AIMAG(this%data_inout(this%bit_rev(j)))  &
                   +wi* REAL (this%data_inout(this%bit_rev(j))), 8 )
        this%data_inout(this%bit_rev(j)) = this%data_inout(this%bit_rev(i)) - tc
        this%data_inout(this%bit_rev(i)) = this%data_inout(this%bit_rev(i)) + tc
      END DO
      wtmp = wr
      wr   = wr*wpr - wi*wpi  +wr
      wi   = wi*wpr + wtmp*wpi+wi
    END DO
    imax = istep
  END DO
  IF ( dir .EQ. FFTHyb_FORWARD ) &
    this%data_inout = this%data_inout*this%fs

END SUBROUTINE FFTHyb_run
!!***

!!****f* ABINIT/m_FFTHyb/FFTHyb_getData
!! NAME
!!  FFTHyb_getData
!!
!! FUNCTION
!!  get result
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=FFT
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

SUBROUTINE FFTHyb_getData(this, bound, array_out, freqs)

!Arguments ------------------------------------
  TYPE(FFTHyb), INTENT(IN) :: this
  INTEGER         ,               INTENT(OUT) :: bound
  DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: array_out
  DOUBLE PRECISION, DIMENSION(:), OPTIONAL, INTENT(OUT) :: freqs
!Local variables ------------------------------
  INTEGER :: i
  INTEGER :: size_out
  DOUBLE PRECISION :: re
  DOUBLE PRECISION :: im

  size_out = SIZE(array_out)
!  IF ( SIZE(array_out) .NE. this%size/2 ) &
!    CALL WARNALL("FFTHyb_getData : size_in != size")
  array_out = 0.d0
  bound = MIN(this%size/2, size_out)
  DO i=1, bound
    re = REAL(this%data_inout(this%bit_rev(i-1)))
    im = AIMAG(this%data_inout(this%bit_rev(i-1)))
    array_out(i) = SQRT(re*re + im*im)
  END DO
  IF ( PRESENT( freqs ) .AND. bound .LE. SIZE(freqs)) THEN
    DO i=1, bound
        freqs(i) = DBLE(i-1)/this%Ts 
    END DO
!  ELSE IF ( PRESENT( freqs ) .AND. bound .GT. SIZE(freqs) ) THEN
!    CALL WARNALL("FFHyb_getData : freqs does is too small")
  END IF
  
END SUBROUTINE FFTHyb_getData
!!***

!!****f* ABINIT/m_FFTHyb/FFTHyb_destroy
!! NAME
!!  FFTHyb_destroy
!!
!! FUNCTION
!!  destroy every thing
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=FFT
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

SUBROUTINE FFTHyb_destroy(this)

!Arguments ------------------------------------
  TYPE(FFTHyb), INTENT(INOUT) :: this

  FREEIF(this%bit_rev)
  FREEIF(this%data_inout)
  this%set = .FALSE.
END SUBROUTINE FFTHyb_destroy
!!***

END MODULE m_FFTHyb
!!***

