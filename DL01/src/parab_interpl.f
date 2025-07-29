      SUBROUTINE PARAB_INTERPL(X,Y,N,T,Z)
      IMPLICIT NONE

c arguments

      INTEGER N
      DOUBLE PRECISION T,Z
      DOUBLE PRECISION X(N),Y(N)

c local variables:

      INTEGER I,J,K,L,M
      DOUBLE PRECISION S
c-----------------------------------------------------------------------
c subroutine parab_interpl
c given:
c     x[] = array of x values
c     y[] = array of y values
c     n   = dimension of x and y
c     t   = desired value f x
c returns
c     z   = interpolated value of y
c
c written by Aigen Li, Princeton University
c history
c 00.11.13 (BTD) cosmetic changes, add comments
c end history
c-----------------------------------------------------------------------

      IF(N.LE.0)THEN
         Z=0
         RETURN
      ELSEIF(N.EQ.1)THEN
         Z=Y(1)
         RETURN
      ELSEIF(N.EQ.2)THEN
         Z=(Y(1)*(T-X(2))-Y(2)*(T-X(1)))/(X(1)-X(2))
         RETURN
      ENDIF

c have eliminated trivial cases. Now deal with n.ge.3:

      IF(T.LE.X(2))THEN
         K=1
         M=3
      ELSEIF(T.GE.X(N-1))THEN
         K=N-2
         M=N
      ELSE
         K=1
         M=N

 10      IF(IABS(K-M).NE.1)THEN
            L=(K+M)/2
            IF(T.LT.X(L))THEN
               M=L
            ELSE
               K=L
            ENDIF
            GOTO 10
         ENDIF
         IF(DABS(T-X(K)).LT.DABS(T-X(M)))THEN
            K=K-1
         ELSE
            M=M+1
         ENDIF
      ENDIF
      Z=0.0
      DO I=K,M
         S=1.0
         DO J=K,M
            IF(J.NE.I)THEN
               S=S*(T-X(J))/(X(I)-X(J))
            ENDIF
         ENDDO
         Z=Z+S*Y(I)
      ENDDO

      RETURN
      END
