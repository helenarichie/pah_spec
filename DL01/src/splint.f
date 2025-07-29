      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      IMPLICIT NONE
C
C Arguments
C
      INTEGER N
      REAL X,Y
      REAL XA(N),Y2A(N),YA(N)
C
C Local variables:
C
      INTEGER K,KHI,KLO
      REAL A,B,H
C***********************************************************************
C Given
C    XA(N) = table of values of independent variable x (monotonic)
C    YA(N) =   "    "  "     "  function y(x) at points XA(n)
C    Y2A(N)= table of second derivatives y"(x) at points XA(n)
C            (can be supplied by routine SPLINE)
C    X     = value of independent variable x
C
C Returns:
C    Y     = estimate for y(X) using cubic spline fit.
C
C This is adapted from Numerical Recipes routine splint.for
C history:
C 03.09.28 (BTD) modified handling of error condition to
C                provide information on array XA
C end history
C***********************************************************************
      KLO=1
      KHI=N
 1000 IF(KHI-KLO.GT.1)THEN
         K=(KHI+KLO)/2
         IF(XA(K).GT.X)THEN
            KHI=K
         ELSE
            KLO=K
         ENDIF
      GOTO 1000
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF(H.EQ.0.)THEN
         WRITE(0,*)'*** Fatal error in subroutine SPLINT: bad XA input'
         WRITE(0,*)'*** KLO,KHI=',KLO,KHI
         WRITE(0,*)'*** XA(KLO)=',XA(KLO),' XA(KHI)=',XA(KHI)
         STOP
      ENDIF
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*
     &   (H**2)/6.
      RETURN
      END
