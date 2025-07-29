      FUNCTION ZBRENT_D(FUNC,X1,X2,TOL)
C Arguments:
      DOUBLE PRECISION FUNC,TOL,X1,X2,ZBRENT_D
      EXTERNAL FUNC
C Parameters:
      DOUBLE PRECISION EPS
      INTEGER ITMAX
      PARAMETER(ITMAX=100,EPS=3.D-8)
      INTEGER ITER
      DOUBLE PRECISION A,B,C,D,E,FA,FB,FC,P,Q,R,S,TOL1,XM
C***********************************************************************
C Purpose: to find zero of function FUNC lying between X1 and X2
C This routine is taken from Numerical Recipes (Press et al.) and
C only slightly modified.
C
C***********************************************************************
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
      IF((FA.GT.0..AND.FB.GT.0.).OR.(FA.LT.0..AND.FB.LT.0.))THEN
         WRITE(0,FMT='(A)')
     &      'Fatal error: root must be bracketed for zbrent'
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
            ZBRENT_D=B
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
         FB=FUNC(B)
      ENDDO
      WRITE(0,FMT='(A)')
     &   'warning: zbrent exceeding maximum iterations'
      ZBRENT_D=B
      RETURN
      END
