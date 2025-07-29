      FUNCTION AGLI_ZBRENT_D(SIZE,NISRF,ISRF_WL,
     &                       DLGLAMBDA,CABS,RATE_HEATING,FUNC,
     &                       X1,X2,TOL,SCRW)

c-----------------       agli_zbrent_d_v2     --------------------
      IMPLICIT NONE

c arguments:

      DOUBLE PRECISION AGLI_ZBRENT_D

      INTEGER NISRF	
      DOUBLE PRECISION DLGLAMBDA,SIZE,RATE_HEATING,X1,X2,TOL
      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF_WL(NISRF),
     &   SCRW(NISRF)

c external functions

      DOUBLE PRECISION FUNC
      EXTERNAL FUNC

c parameters:

      DOUBLE PRECISION EPS
      INTEGER ITMAX
      PARAMETER(ITMAX=100,EPS=3.D-8)

c local variables

      INTEGER ITER
      DOUBLE PRECISION A,B,C,D,E,FA,FB,FC,P,Q,R,S,TOL1,XM

c-----------------------------------------------------------------------
c function AGLI_ZBRENT_D
c given
c     SIZE = grain radius (cm)
c     NISRF = number of wavelengths in tabulated ISRF
c     ISRF_WL[] = wavelengths (cm) for tabulated ISRF
c                 (assumed to be uniform in lg(lambda)
c     DLGLAMBDA = delta(log10(lambda)) for wavelengths ISRF_WL
c     CABS[] = C_abs (cm) at wavelengths ISRF_WL[]
c     RATE_HEATING = heating rate P_heat (erg s-1) for grain in ISRF
c     X1 = lower temperature bracketing solution
c     X2 = upper temperture bracketing solution
c     TOL = fractional error tolerance
c and
c     FUNC = external function to calculate 
c            FUNC=(P_cool-P_heat)/(P_cool+P_heat)
c            where P_cool = radiative cooling rate (erg s-1)
c returns
c     AGLI_ZBRENT_D = temperature (K) at which P_cool = P_heat
c
c Adapted from routine ZBRENT in Numerical Recipes 
c (Press et al.) with only minor modifications.
c Originally adapted by Aigen Li, Princeton University
c History:
c 00.11.12 (BTD) cosmetic changes, added comments
c 00.11.28 (BTD) add SCRW to argument list so scratch space can
c                be passed to function AGLI_FZERO
c 14.02.19 (BTD) v2
c                eliminated irrelevant argument COMPOSITION
c end history
c-----------------------------------------------------------------------
c*** diagnostic
c      write(0,*)'entered function agli_zbrent_d'
c      write(0,*)'X1=',X1,' TOL=',TOL
c***
      A=X1
      B=X2
      FA=FUNC(NISRF,ISRF_WL,DLGLAMBDA,CABS,
     &        RATE_HEATING,A,SCRW)
c*** diagnostic
c      write(0,*)'returned from call A, do call B...'
c***
      FB=FUNC(NISRF,ISRF_WL,DLGLAMBDA,CABS,
     &        RATE_HEATING,B,SCRW)
      IF((FA.GT.0..AND.FB.GT.0.).OR.(FA.LT.0..AND.FB.LT.0.))THEN
         WRITE(0,*)'Fatal error: in agli_zbrent_d: ',
     &             'root must be bracketed for zbrent'
         STOP
      ENDIF
      C=B
      FC=FB
      DO ITER=1,ITMAX
         IF((FB.GT.0..AND.FC.GT.0.).OR.(FB.LT.0..AND.FC.LT.0.))THEN
            C=A
            FC=FA
            D=B-A
            E=D
         ENDIF
         IF(ABS(FC).LT.ABS(FB)) THEN
            A=B
            B=C
            C=A
            FA=FB
            FB=FC
            FC=FA
         ENDIF
         TOL1=2.D0*EPS*ABS(B)+0.5D0*TOL
         XM=.5D0*(C-B)
         IF(ABS(XM).LE.TOL1.OR.FB.EQ.0.)THEN
            AGLI_ZBRENT_D=B
            RETURN
         ENDIF
         IF(ABS(E).GE.TOL1.AND.ABS(FA).GT.ABS(FB))THEN
            S=FB/FA
            IF(A.EQ.C) THEN
               P=2.D0*XM*S
               Q=1.D0-S
            ELSE
               Q=FA/FC
               R=FB/FC
               P=S*(2.D0*XM*Q*(Q-R)-(B-A)*(R-1.D0))
               Q=(Q-1.D0)*(R-1.D0)*(S-1.D0)
            ENDIF
            IF(P.GT.0.)Q=-Q
            P=ABS(P)
            IF(2.D0*P.LT.MIN(3.D0*XM*Q-ABS(TOL1*Q),ABS(E*Q)))THEN
               E=D
               D=P/Q
            ELSE
               D=XM
               E=D
            ENDIF
         ELSE
            D=XM
            E=D
         ENDIF
         A=B
         FA=FB
         IF(ABS(D).GT.TOL1)THEN
            B=B+D
         ELSE
            B=B+SIGN(TOL1,XM)
         ENDIF
         FB=FUNC(NISRF,ISRF_WL,DLGLAMBDA,CABS,
     &           RATE_HEATING,B,SCRW)
      ENDDO
      WRITE(0,*)'Fatal error in agli_zbrent: exceeded max iterations'
      STOP
      END
