!!****m* ABINIT/m_random_zbq
!! NAME
!! m_random_zbq
!!
!! FUNCTION
!!   This library contains routines for generating uniform random numbers
!!   between 0 & 1, using a Marsaglia-Zaman type subtract-with-borrow generator.
!!   This random number generator is non deterministic (it is not "pseudo-random").
!!
!!   These routines have been downloaded from Richard Chandler homepage,
!!   then converted into Fortran90 standard.
!!   Authors:
!!      Richard Chandler (email: richard@stats.ucl.ac.uk)
!!      Paul Northrop (email: paul@stats.ucl.ac.uk)
!!
!! COPYRIGHT
!!   Copyright (C) Richard Chandler/Paul Northrop
!!   Please feel free to use and adapt these routines for your own use. We
!!   do, however, request that if the routines are used in work which is later
!!   presented in public, or published, the source of the code (ie. us!) is
!!   acknowledged. If for no other reason, this will allow other people to test
!!   the routines to see if there's anything horribly wrong with them that will
!!   invalidate your results.
!!
!!   See http://www.ucl.ac.uk/~ucakarc/work/software/randgen.txt
!!   for further information
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


module m_random_zbq

 use m_abicore

 use defs_basis,   only : sp, dp, std_out

 implicit none

 private
!!***

 public :: ZBQLU01
 public :: ZBQLINI

 REAL(DP),SAVE,PRIVATE :: B=4.294967291D9,C=0.0D0
 REAL(DP),SAVE,PRIVATE :: ZBQLIX(43)= &
&     (/8.001441D7,5.5321801D8,&
&       1.69570999D8,2.88589940D8,2.91581871D8,1.03842493D8,&
&       7.9952507D7,3.81202335D8,3.11575334D8,4.02878631D8, &
&       2.49757109D8,1.15192595D8,2.10629619D8,3.99952890D8,&
&       4.12280521D8,1.33873288D8,7.1345525D7,2.23467704D8, &
&       2.82934796D8,9.9756750D7,1.68564303D8,2.86817366D8, &
&       1.14310713D8,3.47045253D8,9.3762426D7 ,1.09670477D8,&
&       3.20029657D8,3.26369301D8,9.441177D6,3.53244738D8,  &
&       2.44771580D8,1.59804337D8,2.07319904D8,3.37342907D8,&
&       3.75423178D8,7.0893571D7 ,4.26059785D8,3.95854390D8,&
&       2.0081010D7,5.9250059D7,1.62176640D8,3.20429173D8,  &
&       2.63576576D8/)


contains

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!*******************************************************************
      REAL(DP) FUNCTION ZBQLU01(DUMMY)
!*******************************************************************
!       Returns a uniform random number between 0 & 1, using
!       a Marsaglia-Zaman type subtract-with-borrow generator.
!       Uses double precision, rather than integer, arithmetic
!       throughout because MZ's integer constants overflow
!       32-bit integer storage (which goes from -2^31 to 2^31).
!       Ideally, we would explicitly truncate all integer
!       quantities at each stage to ensure that the double
!       precision representations do not accumulate approximation
!       error; however, on some machines the use of DNINT to
!       accomplish this is *seriously* slow (run-time increased
!       by a factor of about 3). This double precision version
!       has been tested against an integer implementation that
!       uses long integers (non-standard and, again, slow) -
!       the output was identical up to the 16th decimal place
!       after 10^10 calls, so we're probably OK ...
!
      REAL(DP) :: DUMMY

      INTEGER,SAVE :: CURPOS=1,ID22=22,ID43=43
      REAL(DP) :: X,B2,BINV

      B2 = B
      BINV = 1.0D0/B
 5    X = ZBQLIX(ID22) - ZBQLIX(ID43) - C
      IF (X<0.0D0) THEN
       X = X + B
       C = 1.0D0
      ELSE
       C = 0.0D0
      ENDIF
      ZBQLIX(ID43) = X
!
!     Update array pointers. Do explicit check for bounds of each to
!     avoid expense of modular arithmetic. If one of them is 0 the others
!     won't be
!
      CURPOS = CURPOS - 1
      ID22 = ID22 - 1
      ID43 = ID43 - 1
      IF (CURPOS==0) THEN
       CURPOS=43
      ELSEIF (ID22==0) THEN
       ID22 = 43
      ELSEIF (ID43==0) THEN
       ID43 = 43
      ENDIF
!
!     The integer arithmetic there can yield X=0, which can cause
!     problems in subsequent routines (e.g. ZBQLEXP). The problem
!     is simply that X is discrete whereas U is supposed to
!     be continuous - hence if X is 0, go back and generate another
!     X and return X/B^2 (etc.), which will be uniform on (0,1/B).
      IF (X<BINV) THEN
       B2 = B2*B
       GOTO 5
      ENDIF

      ZBQLU01 = X/B2

      RETURN;X=DUMMY

      END FUNCTION ZBQLU01

!*********************************************************************
      SUBROUTINE ZBQLINI(SEED)
!*********************************************************************
!       To initialize the random number generator - either
!       repeatably or nonrepeatably. Need double precision
!       variables because integer storage can't handle the
!       numbers involved
!*********************************************************************
!	ARGUMENTS
!	=========
!	SEED	(integer, input). User-input number which generates
!		elements of the array ZBQLIX, which is subsequently used
!		in the random number generation algorithm. If SEED=0,
!		the array is seeded using the system clock if the
!		FORTRAN implementation allows it.
!*********************************************************************
      INTEGER :: SEED
!*********************************************************************
!	PARAMETERS
!	==========
!	LFLNO	(integer). Number of lowest file handle to try when
!		opening a temporary file to copy the system clock into.
!		Default is 80 to keep out of the way of any existing
!		open files (although the program keeps searching till
!		it finds an available handle). If this causes problems,
!               (which will only happen if handles 80 through 99 are
!               already in use), decrease the default value.
!*********************************************************************
!      INTEGER,PARAMETER :: LFLNO=80
!*********************************************************************
!	VARIABLES
!	=========
!	SEED	See above
!	ZBQLIX	Seed array for the random number generator. Defined
!		in ZBQLBD01
!	B,C	Used in congruential initialisation of ZBQLIX
!	SS,MM,}	System clock secs, mins, hours and days
!	HH,DD }
!	FILNO	File handle used for temporary file
!	INIT	Indicates whether generator has already been initialised

      INTEGER,SAVE :: INIT
      INTEGER :: SS,MM,HH,DD,I
      REAL(DP) :: TMPVAR1,DSS,DMM,DHH,DDD
!     Variable used by date_and_time function:
      CHARACTER(LEN=8)  :: DATE
      CHARACTER(LEN=10) :: TIME
      CHARACTER(LEN=5)  :: ZONE
      INTEGER :: VALUES(8)

!     Ensure we don't call this more than once in a program
      IF (INIT>=1) THEN
       IF(INIT==1) THEN
         WRITE(std_out,1)
 1       FORMAT(//5X,'****WARNING**** You have called routine ZBQLINI ',&
     &   'more than',/5X,'once. I''m ignoring any subsequent calls.',//)
         INIT = 2
       END IF
       RETURN
      ELSE
       INIT = 1
      ENDIF

!     If SEED = 0, cat the contents of the clock into a file
!     and transform to obtain ZQBLIX(1), then use a congr.
!     algorithm to set remaining elements. Otherwise take
!     specified value of SEED.
!
      IF (SEED==0) THEN
!        Initial coding:
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>	COMMENT OUT FROM HERE IF YOU DON'T HAVE  >>>>>>>>>>>>
!>>>>>>>	'CALL SYSTEM' CAPABILITY ...		 >>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!        CALL SYSTEM(' date +%S%M%H%j > zbql1234.tmp')
! !      Try all file numbers for LFLNO to 999
!        FILNO = LFLNO
!  10    OPEN(FILNO,FILE='zbql1234.tmp',ERR=11)
!        GOTO 12
!  11    FILNO = FILNO + 1
!        IF (FILNO.GT.999) THEN
!         WRITE(std_out,*) 'ZBQLINI: Error on temporary file !'
!         RETURN
!        ENDIF
!        GOTO 10
!  12    READ(FILNO,'(3(I2),I3)') SS,MM,HH,DD
!        CLOSE(FILNO)
!        CALL SYSTEM('rm zbql1234.tmp')

!        New F90 coding:
         CALL DATE_AND_TIME(DATE,TIME,ZONE,VALUES)
         HH=VALUES(4)
         MM=VALUES(5)
         SS=VALUES(6)
         DD=(VALUES(2)-1)*30+VALUES(3) ! Approx. for the day index

         DSS = DINT((DBLE(SS)/6.0D1) * B)
         DMM = DINT((DBLE(MM)/6.0D1) * B)
         DHH = DINT((DBLE(HH)/2.4D1) * B)
         DDD = DINT((DBLE(DD)/3.65D2) * B)
         TMPVAR1 = DMOD(DSS+DMM+DHH+DDD,B)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<	... TO HERE (END OF COMMENTING OUT FOR 	  <<<<<<<<<<<
!<<<<<<<<	USERS WITHOUT 'CALL SYSTEM' CAPABILITY	  <<<<<<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ELSE
         TMPVAR1 = DMOD(DBLE(SEED),B)
      ENDIF

      ZBQLIX(1) = TMPVAR1
      DO I = 2,43
       TMPVAR1 = ZBQLIX(I-1)*3.0269D4
       TMPVAR1 = DMOD(TMPVAR1,B)
       ZBQLIX(I) = TMPVAR1
      END DO

      END SUBROUTINE ZBQLINI

end module m_random_zbq
!!***
