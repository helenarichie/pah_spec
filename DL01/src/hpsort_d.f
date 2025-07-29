      SUBROUTINE HPSORT(N,RA)
      IMPLICIT NONE
C
C Arguments:
C
      INTEGER N
      DOUBLE PRECISION RA(N)
C
C Local variables:
C
      INTEGER I,IR,J,L
      DOUBLE PRECISION RRA
C-----------------------------------------------------------------------
C Subroutine HPSORT
C Given:
C     N = length of array RA
C     RA(1-N) = array values
C Returns:
C     RA(1-N) = array values sorted into ascending numerical order
C
C HPSORT, adapted from Press et al (1992), Numerical Recipes,
C uses the heap sort algorithm
C requires N*log_2(N) operations to complete sorting in worst case
C 
C History:
C 01.06.05 (BTD) added comments
C end history
C-----------------------------------------------------------------------
      IF(N.LT.2)RETURN
      L=N/2+1
      IR=N
10    CONTINUE
      IF(L.GT.1)THEN
         L=L-1
         RRA=RA(L)
      ELSE
         RRA=RA(IR)
         RA(IR)=RA(1)
         IR=IR-1
         IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
         ENDIF
      ENDIF
      I=L
      J=L+L
20    IF(J.LE.IR)THEN
         IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
         ENDIF
         IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
         ELSE
            J=IR+1
         ENDIF
      GOTO 20
      ENDIF
      RA(I)=RRA
      GOTO 10
      END

