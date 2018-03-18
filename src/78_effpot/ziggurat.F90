! Marsaglia & Tsang generator for random normals & random exponentials.
! Translated from C by Alan Miller (amiller@bigpond.net.au)

! Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
! random variables', J. Statist. Software, v5(8).

! This is an electronic journal which can be downloaded from:
! http://www.jstatsoft.org/v05/i08

! N.B. It is assumed that all integers are 32-bit.
! N.B. The value of M2 has been halved to compensate for the lack of
!      unsigned integers in Fortran.

! Latest version - 1 January 2001
MODULE Ziggurat
   IMPLICIT NONE

   PRIVATE

   INTEGER,  PARAMETER  ::  DP=SELECTED_REAL_KIND( 12, 60 )
   REAL(DP), PARAMETER  ::  m1=2147483648.0_DP,   m2=2147483648.0_DP,      &
                            half=0.5_DP
   REAL(DP)             ::  dn=3.442619855899_DP, tn=3.442619855899_DP,    &
                            vn=0.00991256303526217_DP,                     &
                            q,                    de=7.697117470131487_DP, &
                            te=7.697117470131487_DP,                       &
                            ve=0.003949659822581572_DP
   INTEGER,  SAVE       ::  iz, jz, jsr=123456789, kn(0:127),              &
                            ke(0:255), hz
   REAL(DP), SAVE       ::  wn(0:127), fn(0:127), we(0:255), fe(0:255)
   LOGICAL,  SAVE       ::  initialized=.FALSE.

   PUBLIC  :: zigset, shr3, uni, rnor, rexp


CONTAINS


SUBROUTINE zigset( jsrseed )

   INTEGER, INTENT(IN)  :: jsrseed

   INTEGER  :: i

   !  Set the seed
   jsr = jsrseed

   !  Tables for RNOR
   q = vn*EXP(half*dn*dn)
   kn(0) = (dn/q)*m1
   kn(1) = 0
   wn(0) = q/m1
   wn(127) = dn/m1
   fn(0) = 1.0_DP
   fn(127) = EXP( -half*dn*dn )
   DO  i = 126, 1, -1
      dn = SQRT( -2.0_DP * LOG( vn/dn + EXP( -half*dn*dn ) ) )
      kn(i+1) = (dn/tn)*m1
      tn = dn
      fn(i) = EXP(-half*dn*dn)
      wn(i) = dn/m1
   END DO

   !  Tables for REXP
   q = ve*EXP( de )
   ke(0) = (de/q)*m2
   ke(1) = 0
   we(0) = q/m2
   we(255) = de/m2
   fe(0) = 1.0_DP
   fe(255) = EXP( -de )
   DO  i = 254, 1, -1
      de = -LOG( ve/de + EXP( -de ) )
      ke(i+1) = m2 * (de/te)
      te = de
      fe(i) = EXP( -de )
      we(i) = de/m2
   END DO
   initialized = .TRUE.
   RETURN
END SUBROUTINE zigset



!  Generate random 32-bit integers
FUNCTION shr3( ) RESULT( ival )
   INTEGER  ::  ival

   jz = jsr
   jsr = IEOR( jsr, ISHFT( jsr,  13 ) )
   jsr = IEOR( jsr, ISHFT( jsr, -17 ) )
   jsr = IEOR( jsr, ISHFT( jsr,   5 ) )
   ival = jz + jsr
   RETURN
END FUNCTION shr3



!  Generate uniformly distributed random numbers
FUNCTION uni( ) RESULT( fn_val )
   REAL(DP)  ::  fn_val

   fn_val = half + 0.2328306e-9_DP * shr3( )
   RETURN
END FUNCTION uni



!  Generate random normals
FUNCTION rnor( ) RESULT( fn_val )
   REAL(DP)             ::  fn_val

   REAL(DP), PARAMETER  ::  r = 3.442620_DP
   REAL(DP)             ::  x, y

   IF( .NOT. initialized ) CALL zigset( jsr )
   hz = shr3( )
   iz = IAND( hz, 127 )
   IF( ABS( hz ) < kn(iz) ) THEN
      fn_val = hz * wn(iz)
   ELSE
      DO
         IF( iz == 0 ) THEN
            DO
               x = -0.2904764_DP* LOG( uni( ) )
               y = -LOG( uni( ) )
               IF( y+y >= x*x ) EXIT
            END DO
            fn_val = r+x
            IF( hz <= 0 ) fn_val = -fn_val
            RETURN
         END IF
         x = hz * wn(iz)
         IF( fn(iz) + uni( )*(fn(iz-1)-fn(iz)) < EXP(-half*x*x) ) THEN
            fn_val = x
            RETURN
         END IF
         hz = shr3( )
         iz = IAND( hz, 127 )
         IF( ABS( hz ) < kn(iz) ) THEN
            fn_val = hz * wn(iz)
            RETURN
         END IF
      END DO
   END IF
   RETURN
END FUNCTION rnor



!  Generate random exponentials
FUNCTION rexp( ) RESULT( fn_val )
   REAL(DP)  ::  fn_val

   REAL(DP)  ::  x

   IF( .NOT. initialized ) CALL Zigset( jsr )
   jz = shr3( )
   iz = IAND( jz, 255 )
   IF( ABS( jz ) < ke(iz) ) THEN
      fn_val = ABS(jz) * we(iz)
      RETURN
   END IF
   DO
      IF( iz == 0 ) THEN
         fn_val = 7.69711 - LOG( uni( ) )
         RETURN
      END IF
      x = ABS( jz ) * we(iz)
      IF( fe(iz) + uni( )*(fe(iz-1) - fe(iz)) < EXP( -x ) ) THEN
         fn_val = x
         RETURN
      END IF
      jz = shr3( )
      iz = IAND( jz, 255 )
      IF( ABS( jz ) < ke(iz) ) THEN
         fn_val = ABS( jz ) * we(iz)
         RETURN
      END IF
   END DO
   RETURN
END FUNCTION rexp

END MODULE ziggurat



! PROGRAM test_ziggurat
!    USE Ziggurat
!    IMPLICIT NONE

!    INTEGER, PARAMETER  ::  DP = SELECTED_REAL_KIND(12, 60)
!    REAL                ::  t1, t2
!    REAL(DP)            ::  x, chisq, expctd
!    INTEGER             ::  i, ten_million=10000000, freq(0:999), j

!    !  Set the random number seed
!    WRITE(*, '(a)', ADVANCE='NO') ' Enter a non-zero integer: '
!    READ(*, *) i

!    CALL zigset( i )

!    !  First time the generation of uniform, normal & exponential generators
!    CALL CPU_TIME( t1 )
!    DO  i = 1, ten_million
!       x = uni( )
!    END DO
!    CALL CPU_TIME( t2 )
!    WRITE(*, '(a, i10, a, f9.2, a)') ' Time to generate ', ten_million,  &
!                                     ' uniform random nos. = ', t2-t1, 'sec.'
!    CALL CPU_TIME( t1 )
!    DO  i = 1, ten_million
!       x = rnor( )
!    END DO
!    CALL CPU_TIME( t2 )
!    WRITE(*, '(a, i10, a, f9.2, a)') ' Time to generate ', ten_million,     &
!                                     ' normal random nos. = ', t2-t1, 'sec.'

!    CALL CPU_TIME( t1 )
!    DO  i = 1, ten_million
!       x = rexp( )
!    END DO
!    CALL CPU_TIME( t2 )
!    WRITE(*, '(a, i10, a, f9.2, a)') ' Time to generate ', ten_million,     &
!                                     ' exponential random nos. = ', t2-t1, 'sec.'
!    WRITE(*, *)

!    !  Now look at the distributions generated
!    WRITE(*, *) 'Std. devn. of chi-squared with 999 d.of.f. = 44.70'
!    WRITE(*, *) 'Values of chi-squared below should be within about 90. of 999.'

!    freq = 0
!    DO  i = 1, ten_million
!       j = 1000 * uni( )
!       freq(j) = freq(j) +1
!    END DO
!    chisq = 0.0_DP
!    expctd = ten_million / 1000
!    DO  j = 0, 999
!       chisq = chisq + (freq(j) - expctd)**2 / expctd
!    END DO
!    WRITE(*, '(a, f9.2, a)') ' Uniform distribution, chi-squared = ', chisq,   &
!                             ' with 999 deg. of freedom'
!    freq = 0
!    DO  i = 1, ten_million
!       j = 1000 * alnorm( rnor( ), .FALSE. )
!       freq(j) = freq(j) +1
!    END DO
!    chisq = 0.0_DP
!    expctd = ten_million / 1000
!    DO  j = 0, 999
!       chisq = chisq + (freq(j) - expctd)**2 / expctd
!    END DO
!    WRITE(*, '(a, f9.2, a)') ' Normal distribution, chi-squared = ', chisq,  &
!                             ' with 999 deg. of freedom'
!    freq = 0
!    DO  i = 1, ten_million
!       j = 1000 * EXP( -rexp( ) )
!       IF( j > 999  .OR.  j < 0 ) THEN
!          WRITE(*, '(a, 2i10)') ' i, j = ', i, j
!       ELSE
!          freq(j) = freq(j) +1
!       END IF
!    END DO
!    chisq = 0.0_DP
!    expctd = ten_million / 1000
!    DO  j = 0, 999
!       chisq = chisq + (freq(j) - expctd)**2 / expctd
!    END DO
!    WRITE(*, '(a, f9.2, a)') ' Exponential distribution, chi-squared = ',      &
!                             chisq, ' with 999 deg. of freedom'
!    STOP



! CONTAINS


! !  Algorithm AS66 Applied Statistics (1973) vol.22, no.3

! !  Evaluates the tail area of the standardised normal curve
! !  from x to infinity if upper is .true. or
! !  from minus infinity to x if upper is .false.

! ! ELF90-compatible version by Alan Miller
! ! Latest revision - 29 November 2001

! FUNCTION alnorm( x, upper ) RESULT( fn_val )
!    REAL(DP), INTENT(IN)   ::  x
!    LOGICAL,   INTENT(IN)  ::  upper
!    REAL(DP)               ::  fn_val

!    !  Local variables
!    REAL(DP), PARAMETER   ::  zero=0.0_DP, one=1.0_DP, half=0.5_DP, con=1.28_DP
!    REAL(DP)              ::  z, y
!    LOGICAL               ::  up

!    !  Machine dependent constants
!    REAL(DP), PARAMETER  ::  ltone = 7.0_DP, utzero = 18.66_DP
!    REAL(DP), PARAMETER  ::  p = 0.398942280444_DP, q = 0.39990348504_DP,   &
!                             r = 0.398942280385_DP, a1 = 5.75885480458_DP,  &
!                             a2 = 2.62433121679_DP, a3 = 5.92885724438_DP,  &
!                             b1 = -29.8213557807_DP, b2 = 48.6959930692_DP, &
!                             c1 = -3.8052E-8_DP, c2 = 3.98064794E-4_DP,     &
!                             c3 = -0.151679116635_DP, c4 = 4.8385912808_DP, &
!                             c5 = 0.742380924027_DP, c6 = 3.99019417011_DP, &
!                             d1 = 1.00000615302_DP, d2 = 1.98615381364_DP,  &
!                             d3 = 5.29330324926_DP, d4 = -15.1508972451_DP, &
!                             d5 = 30.789933034_DP

!    up = upper
!    z = x
!    IF( z < zero ) THEN
!       up = .NOT. up
!       z = -z
!    END IF
!    IF( z <= ltone  .OR.  (up  .AND.  z <= utzero) ) THEN
!       y = half*z*z
!       IF( z > con ) THEN
!          fn_val = r*EXP( -y )/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))))
!       ELSE
!          fn_val = half - z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))))
!       END IF
!    ELSE
!       fn_val = zero
!    END IF

!    IF( .NOT. up ) fn_val = one - fn_val
!    RETURN
! END FUNCTION alnorm

! END PROGRAM test_ziggurat

