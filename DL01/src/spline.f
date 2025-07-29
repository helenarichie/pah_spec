      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      IMPLICIT NONE
C
C Arguments
C
      INTEGER N
      REAL YP1,YPN
      REAL X(N),Y(N),Y2(N)
C
C Local variables
C
      INTEGER NMAX
      PARAMETER(NMAX=1000)
      INTEGER I,K
      REAL P,QN,SIG,UN
      REAL U(NMAX)
C***********************************************************************
C Given:
C    N   = number of data points
C    X(N)=table of data points
C    Y(N)=  "    "  "     "
C    YP1 = first derivative dY/dX at X(1)
C    YPN = first derivative dY/dX at X(N)
C
C Returns:
C    Y2(N)=table of second derivatives of the interpolating function
C          at points X(1)-X(N)
C
C If YP1 and/or YPN are equal to 1.e30 or larger, the routine will
C set the corresponding boundary condition for a "natural" spline,
C with zero second derivative on that boundary.
C
C This is routine spline.for from Numerical Recipes.
C 
C***********************************************************************

c sanity check
      IF(N.GT.NMAX)THEN
         WRITE(0,*)'Need to recompile spline.f with NMAX.ge.',N
         STOP
      ENDIF
C
C Determine whether to use "natural" boundary condition or
C specified first derivatives at boundaries.
C
      IF(YP1.GT.0.99E30)THEN
         Y2(1)=0.
         U(1)=0.
      ELSE
         Y2(1)=-0.5
         U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO I=2,N-1
         SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
         P=SIG*Y2(I-1)+2.
         Y2(I)=(SIG-1.)/P
         U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))/
     &        (X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      ENDDO
c*** diagnostic
c      write(0,*)'in spline, N=',N
c      write(0,*)'U(I=1,N-1)=',(U(K),K=1,N-1)
c      write(0,*)'U(52)=',U(52)
c      write(0,*)'U(53)=',U(53)
c      write(0,*)'X(50-54)=',(X(K),K=50,54)
c      write(0,*)'Y(50-54)=',(Y(K),K=50,54)
c***
C
C Determine whether to use "natural" boundary condition or
C specified first derivatives at boundaries.
C
      IF(YPN.GT.0.99E30)THEN
         QN=0.
         UN=0.
      ELSE
         QN=0.5
         UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
C
C Backsubstitution loop of tridiagonal algorithm:
C
      DO K=N-1,1,-1
         Y2(K)=Y2(K)*Y2(K+1)+U(K)
      ENDDO
      RETURN
      END
