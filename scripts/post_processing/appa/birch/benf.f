      PROGRAM EXBIRFT
C     Fits energy versus volume data to the form
C     E(V) = Sum [ a(n) v^(-2n/3) , {n,0,N}]
C     and then calculates E(Vo), Vo, Ko, and Ko'

      PARAMETER (MAXN=50,MAXM=10,MAXIN=200)
      IMPLICIT DOUBLE PRECISION (A-H,K,O-Z)
      DIMENSION XIN(MAXIN)
      LOGICAL CONA,DISB,ENGR,MINFOUND,FIRST
      DIMENSION V(MAXN),E(MAXN),PR(MAXN),P(0:MAXM,MAXN)
      DIMENSION CONV(4),B(0:MAXM)
      common /ear/ a(0:maxm),earray(0:3),m,n
      external efun,defun
      CHARACTER*50 DATAFL,OUTFL,LAB*80
      CHARACTER*1 BLANK,ANS,YESU,YESL
      REAL*4 gamma0
C! Conversion factor:  V=CONV*a*a*a
C!    1:sc, 2:fcc, 3:bcc, 4:hex
      DATA CONV/1.d0,2.5D-1,5D-1,6.876241706D-1/
C! 1 au in Angstroms
      DATA AUOAN/5.2917706D-1/
C! 1 Ry in eV
      DATA RYOEV/1.3605804D1/

      DATA BLANK/' '/,YESU/'Y'/,YESL/'y'/
      DATA THIRD/3.33333333333333333D-1/
      DATA TWO3/6.666666666666666667D-1/

      AU3=AUOAN*AUOAN*AUOAN

      OPEN(UNIT=21,FILE='fit.res')
      OPEN(UNIT=12,FILE='fit.dat')
      OPEN(UNIT=11,FILE='fit.inp',STATUS='OLD')
 3000 continue

      FIRST=.TRUE.
      READ (11,*)NCOL
 2000 CONTINUE
      READ(11, *)IX
      IF(IX.EQ.0) THEN
       CLOSE(21)
       CLOSE(12)
       CLOSE(11)
       STOP
      ENDIF

      READ (11,*)IY
      N=0
      IF(FIRST)THEN
       FIRST=.FALSE.
       READ(11,*) IWHICH
       CONA=IWHICH.EQ.2
       IF(.NOT.CONA) IWHICH=1
       READ(11,*) LAT
       READ(11,*) IAUNIT
       DISB=IAUNIT.EQ.2
       IF(.NOT.DISB) IAUNIT=1
       IF(DISB) THEN
         AUNIT=1D0
       ELSE
         AUNIT=AUOAN
       END IF
       VUNIT=AUNIT*AUNIT*AUNIT
      READ(11,*) IEUNIT
      IF(IEUNIT.EQ.1) THEN
         EUNIT=1D0
       ELSE IF(IEUNIT.EQ.2) THEN
          EUNIT=RYOEV
       ELSE IF(IEUNIT.EQ.3) THEN
          EUNIT=0.5D0
       ELSE IF(IEUNIT.EQ.4)THEN
          EUNIT=1312.0125D0
       END IF
      ENDIF
 5000 CONTINUE
      READ(11,*) CMULTV
      READ(11,*) CMULTE
C
C     To find the minimum energy, we'll need an estimated starting value.
C       Since we've got the energies on file, we might as well use the
C       volume of the lowest energy as our estimated Vo.
      EMIN=1D10
      VMIN=0D0
C     Read data in the form X,E (x=volume or lattice constant, E=energy)
C     where Angstroms and eV or atomic units (Bohrs and Rydbergs) are
C     used depending upont the setting of IUNIT
      N=0
      READ(11,*)npoints
      do np=1,npoints
200   READ(11,*)(XIN(I),I=1,NCOL)
      N=N+1
      X=XIN(IX)
      EI=XIN(IY)
      EE=XIN(IY)/CMULTE
      IF(CONA) THEN
        VA=CONV(LAT)*X*X*X*cmultV
      ELSE
        VA=X*cmultV
      END IF
      V(N)=VA/VUNIT
C! Volume in au**3
      E(N)=EE/EUNIT
C! Energy in Ryd
      WRITE(21,256) N,X,VA,EI,EE
256   FORMAT(1X,I5,4(G15.9,2x))
      IF(E(N).LT.EMIN) THEN
        EMIN=E(N)
        VMIN=V(N)
        nsave=n
      END IF
      enddo
c      GO TO 200
300   CONTINUE
C     XSTART is the starting value for the minimum search:
      XSTART=VMIN**(-TWO3)
      WRITE(21,305) EMIN,VMIN,XSTART
305   FORMAT(/' Minimum energy in file = ',F15.5,' at V = ',
     1        2F15.7/)
C
C     Set up the fit.  Note that the functions to be fitted are
C       v^(-2n/3),n=0,1,2,...MAXM
4000  continue
      READ(11,*) M
      IF(M.eq.0) GO TO 5000
      IF(M.lt.0) GO TO 3000
      IF(M.GT.MAXM) GO TO 4000
      write(21,29)M
  29  Format(/'Order of v^(-2/3) polynomial: ',i1,/)
      M1=M+1
      DO 400 I=1,N

C       Establish the basis functions:
        P(0,I)=1D0
        X=V(I)**(-TWO3)
        DO 400 J=1,M
400       P(J,I)=X**J
      CALL LSTSQR(N,M1,E,A,P)
      do I=0,M
      B(I)=A(I)*1312.0125D0
      enddo
      WRITE(21,415) (I,A(I),I=0,M)
      WRITE(21,416) (I,B(I),I=0,M)
415   FORMAT(/' Fitting coefficients Ryd   :'/(1X,I5,1PE17.9))
416   FORMAT(/' Fitting coefficients kJ/mol:'/(1X,I5,1PE17.9))
      write(21,417)
417   FORMAT(/)
C
C     Now for the error checking:
C
      ERMS=0D0
      EMAX=0D0
      DO 600 I=1,N
        XO=V(I)**(-TWO3)
        CALL PLYEVL(M,0D0,A,XO,0,EARRAY)
        ECK=EARRAY(0)
        ERR=ECK-E(I)
        IF(ABS(ERR).GT.ABS(EMAX)) EMAX=ERR
        ERMS=ERMS+ERR*ERR
        WRITE(21,605) I,V(I),E(I),ECK,ERR
605     FORMAT(1X,I5,F12.5,3F15.9)
600   continue
      ERMS=SQRT(ERMS/N)
C     Convert the errors to eV
      ERMSEV=ERMS*RYOEV
      EMAXEV=EMAX*RYOEV
C     Now we must find the equilibrium volume, VO, if we can.
      if(nsave.lt.n)then
      vp1=v(nsave+1)
      vm1=v(nsave-1)
      else if(nsave.eq.1)then
      vp1=v(nsave+1)
      vm1=vmin*0.9
      else
      vp1=vmin*1.4
      vm1=v(nsave-1)
      endif
      fO= dbrent(vm1,vmin,vp1,efun,defun,1d-12,vO)
      XO=Vo**(-TWO3)
      CALL PLYEVL(M,0D0,A,Xo,3,EARRAY)
      eo=earray(0)

      p0=TWO3*Xo*EARRAY(1)/vo
      WRITE(21,715) P0
715   FORMAT(' "Equilibrium" pressure = ',1PE15.6)
      KO=(1D1*EARRAY(1)+4D0*XO*EARRAY(2))*XO/(9D0*VO)
      KOP=(2.5D1*EARRAY(1)+4D0*XO*(6D0*EARRAY(2)+XO*EARRAY(3)))/
     1    (3D0*(5D0*EARRAY(1)+2D0*XO*EARRAY(2)))
      KPP=0D0
      gamma0=-1D0+((2D1*EARRAY(1)+1.5D1*EARRAY(2)*XO+2D0*EARRAY(3)
     .*XO*XO)/(5D0*EARRAY(1)+2D0*EARRAY(2)*XO))*(1D0/3D0)
      DO 720 I=0,M
      XI=dfloat(I)
      VX=VO**(-(3D0+TWO3*XI))
      KPP=KPP+TWO3*XI*VX*A(I)*(2D0+TWO3*XI)*(1+TWO3*XI)*(1+TWO3*XI)
 720  CONTINUE
      KPP=KPP*VO*VO/(KO*KO)-(KOP/KO)*(1+KOP)
      ALAT=(VO/CONV(LAT))**3.33333333333333333D-1
      EOkJ=EO*1312.7341D0
      WRITE(21,835) VO,ALAT,EO,EOkJ,KO,KOP,KPP,gamma0,ERMS,EMAX
835   FORMAT(/' Equilibrium parameters for the Birch-Murnaghan',
     1          ' equation:'
     2/'    Vo = ',F15.5,' Bohr**3  a =',F15.5,' Bohrs',
     3/'    Eo = ',F22.12,' Rydbergs'
     8/'    EokJ=',f15.5,' kJ/mol'
     4/'    Ko = ',F15.5,' Rydbergs/Bohr**3'
     5/'    Ko''= ',F15.5,' Ko" = ',F15.5,
     6/'    gruneisen const=',f12.7,
     7/'    RMS error in energy fit = ',F15.5,' Rydbergs'
     8/'    Maximum error in energy fit = ',F15.5,' Rydbergs'/)

      write(12,*)EOkJ
      write(12,*)VO
      write(12,*)B(0)
      write(12,*)B(1)
      write(12,*)B(2)
      if(M.eq.3)then
      write(12,*)B(3)
      endif

C     Now convert to standard units:
      EO=EO*RYOEV
C! In eV
      KO=1.4710756D4*KO
      KPP=KPP/1.4710756D4
C! In GPa
      VOa=VO*AU3
C! In Angstroms **3
      ALAT=ALAT*AUOAN
      WRITE(21,855) VOa,ALAT,EO,KO,KOP,KPP,ERMSEV,EMAXEV
855   FORMAT(/' Equilibrium parameters for the Birch-Murnaghan',
     1          ' equation:'
     2/'    Vo = ',F15.5,' Angstroms**3   a =',F15.5,' Angstroms',
     3/'    Eo = ',F22.12,' eV'
     4/'    Ko = ',F15.5,' GPa'
     5/'    Ko''= ',F15.5,' Ko" = ',F15.5,
     6/'    RMS error in energy fit = ',F15.5,' eV'
     7/'    Maximum error in energy fit = ',F15.5,' eV'/)
255   FORMAT(1X,I5,3G15.7)

C!    WRITE PRESSURE
      IF (M.EQ.3) THEN
       WRITE(21,'(/,"     Volume (Ang^3) vs Pressure (GPa)")')
       DO 615 I=1,N       
        X=(1.D0/V(I))**THIRD
        PR(I)=(TWO3*A(1)*X**5+2.d0*TWO3*A(2)*X**7+2.D0*A(3)*X**9)
     c       *294210.d0/20.d0
        WRITE(21,'(2F22.12)') PR(I),V(I)*AU3
615    CONTINUE
      ENDIF
      
      GO TO 4000
 1000  STOP
       END
C_______________________________________________________________
C
      SUBROUTINE LSTSQR(NN,MM,F,A,P)
C     FITS THE FUNCTION F(X), DEFINED AT THE N POINTS X(I)
C     ( F(I)=F(X(I)) ), TO THE M FUNCTIONS P(I,X)
C     ( P(I,J)=P(I,X(J)) ), USING A LINEARIZED LEAST SQUARES
C     FITTING ROUTINE.  THE COEFFICIENT OF P(I,X) WILL BE
C     RETURNED AS A(I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXN=50,MAXM=11)
      DIMENSION X(MAXN),F(MAXN),A(MAXM),P(MAXM,MAXN)
      DIMENSION G(MAXM,MAXN),AVG2(MAXM),B(MAXM),
     1  C(MAXM,MAXM),SGP(MAXM),D(MAXM,MAXM)
      M=MM
      N=NN
      RN=1D0/DBLE(N)
C     THE G'S LINEAR COMBINATIONS OF THE P'S, CHOSEN BY THE
C     GRAM-SCHMIDT ORTHOGONALIZATION PROCEDURE SO THAT
C     <G(M)*G(N)>=0 IF M.NE.N, WHERE
C     <H>=[H(X(1))+H(X(2))+H(X(3))+...+H(X(N))]/N
C     FOR ANY FUNCTION H(X)
C     CALCULATE THE ITH FUNCTION
      DO 40 I=1,M
         IM=I-1
         IP=I+1
         SG2=0D0
         DO 10 K=IP,M
10         SGP(K)=0D0
C        AT THE JTH POINT
         DO 30 J=1,N
           SUM=0D0
           DO 20 K=1,IM
20           SUM=SUM+C(I,K)*G(K,J)
           G(I,J)=P(I,J)+SUM
           SG2=SG2+G(I,J)*G(I,J)
           DO 30 K=IP,M
30           SGP(K)=SGP(K)+P(K,J)*G(I,J)
C        AVG2(I)=<G(I)*G(I)>
         AVG2(I)=RN*SG2
C        C(K,I)= -<P(K)*G(I)>/<G(I)*G(I)>
         DO 40 K=IP,M
40         C(K,I)=-SGP(K)/SG2
C     SINCE <G(I)*G(J)>=0 FOR I.NE.J, IT'S TRIVIAL TO FIND
C     THE COEFFICIENTS FOR A LEAST SQUARES FIT OF F TO THE G'S
      DO 60 I=1,M
         SUM=0D0
         DO 50 J=1,N
50         SUM=SUM+G(I,J)*F(J)
60       B(I)=RN*SUM/AVG2(I)
C     TO CONVERT THE B'S INTO A'S, WE FIRST NEED TO FIND THE
C     COEFFICIENTS FOR EXPANDING THE G'S IN TERMS OF THE P'S
      DO 80 I=1,M
         D(I,I)=1D0
         IM=I-1
         DO 80 K=1,IM
           SUM=0D0
           DO 70 L=K,IM
70           SUM=SUM+C(I,L)*D(L,K)
80         D(I,K)=SUM
C     FINALLY, WE CAN CHANGE THE B'S INTO A'S
      DO 100 I=1,M
         SUM=0D0
         DO 90 J=I,M
90         SUM=SUM+B(J)*D(J,I)
100      A(I)=SUM
      RETURN
      END
C_______________________________________________________________________
C
      SUBROUTINE  P L Y E V L (M,X0,A,X,N,P)
C     Evaluates the polynomial given by
C
C       p(x) = Sum [ a(i) (x-x0)^i ,{i,0,m}]
C
C       and its first N derivatives
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     Note the dummy indexing on A and P:
      DIMENSION A(0:1),P(0:1)
      Y=X-X0
C     Zeroth order term (in Y)
      IPROD=1
      DO 100 J=0,N
        P(J)=IPROD*A(M)
100     IPROD=IPROD*(M-J)
      DO 200 I=M-1,0,-1
        IPROD=1
        DO 200 J=0,N
          IF(IPROD.GT.0D0) P(J)=P(J)*Y+IPROD*A(I)
200       IPROD=IPROD*(I-J)
      RETURN
      END
      real*8 function efun(v)
      implicit real*8 (a-h,o-z)
      PARAMETER (MAXN=50,MAXM=10,MAXIN=200)
      common /ear/ a(0:maxm),earray(0:3),m,n
      data two3 /0.6666666666666d0/
      X=V**(-TWO3)
C
C       Call the polynomial to evaluate the energy and pressure:
C
       CALL PLYEVL(M,0D0,A,X,1,EARRAY)
       efun=EARRAY(0)
       return
       end
      real*8 function defun(v)
      PARAMETER (MAXN=50,MAXM=10,MAXIN=200)
      implicit real*8 (a-h,o-z)
      common /ear/ a(0:maxm),earray(0:3),m,n
      data two3 /0.6666666666666d0/
      X=V**(-TWO3)
      defun=-TWO3*X*EARRAY(1)/v
      return
      end
      FUNCTION DBRENT(AX,BX,CX,F,DF,TOL,XMIN)
      implicit real*8 (a-h,o-z)
      PARAMETER (ITMAX=1000,ZEPS=1.0E-12)
      LOGICAL OK1,OK2
      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      X=V
      E=0.
      FX=F(X)
      FV=FX
      FW=FX
      DX=DF(X)
c      print *,x,fx,dx
      DV=DX
      DW=DX
      DO 11 ITER=1,ITMAX
        XM=0.5*(A+B)
        TOL1=TOL*ABS(X)+ZEPS
        TOL2=2.*TOL1
        IF(ABS(X-XM).LE.(TOL2-.5*(B-A))) GOTO 3
        IF(ABS(E).GT.TOL1) THEN
          D1=2.*(B-A)
          D2=D1
          IF(DW.NE.DX) D1=(W-X)*DX/(DX-DW)
          IF(DV.NE.DX) D2=(V-X)*DX/(DX-DV)
          U1=X+D1
          U2=X+D2
          OK1=((A-U1)*(U1-B).GT.0.).AND.(DX*D1.LE.0.)
          OK2=((A-U2)*(U2-B).GT.0.).AND.(DX*D2.LE.0.)
          OLDE=E
          E=D
          IF(.NOT.(OK1.OR.OK2))THEN
            GO TO 1
          ELSE IF (OK1.AND.OK2)THEN
            IF(ABS(D1).LT.ABS(D2))THEN
              D=D1
            ELSE
              D=D2
            ENDIF
          ELSE IF (OK1)THEN
            D=D1
          ELSE
            D=D2
          ENDIF
          IF(ABS(D).GT.ABS(0.5*OLDE))GO TO 1
          U=X+D
          IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
          GOTO 2
        ENDIF
1       IF(DX.GE.0.) THEN
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=0.5*E
2       IF(ABS(D).GE.TOL1) THEN
          U=X+D
          FU=F(U)
        ELSE
          U=X+SIGN(TOL1,D)
          FU=F(U)
          IF(FU.GT.FX)GO TO 3
        ENDIF
        DU=DF(U)
c      print *,x,fx,dx
        IF(FU.LE.FX) THEN
          IF(U.GE.X) THEN
            A=X
          ELSE
            B=X
          ENDIF
          V=W
          FV=FW
          DV=DW
          W=X
          FW=FX
          DW=DX
          X=U
          FX=FU
          DX=DU
        ELSE
          IF(U.LT.X) THEN
            A=U
          ELSE
            B=U
          ENDIF
          IF(FU.LE.FW .OR. W.EQ.X) THEN
            V=W
            FV=FW
            DV=DW
            W=U
            FW=FU
            DW=DU
          ELSE IF(FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
            V=U
            FV=FU
            DV=DU
          ENDIF
        ENDIF
11    CONTINUE
      print *,'DBRENT exceeded maximum iterations.'
      print *,'ctrl+C to stop; any other key to continue...'
      read *
c     PAUSE 'DBRENT exceeded maximum iterations.'
3     XMIN=X
      DBRENT=FX
      RETURN
      END

