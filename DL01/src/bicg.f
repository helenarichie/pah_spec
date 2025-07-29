      SUBROUTINE LINBCG(NMAX,SA,IJA,N,B,X,ITOL,TOL,ITMAX,ITER,ERR)
      IMPLICIT NONE

      INTEGER NMAX0
      DOUBLE PRECISION EPS
C	PARAMETER (NMAX0=1024,EPS=1.0D-14)
      PARAMETER (NMAX0=3024,EPS=1.0D-14)
	
      INTEGER NMAX,N,ITOL,ITMAX,ITER,IJA(NMAX)
      DOUBLE PRECISION SA(NMAX),B(N),X(N),TOL,ERR

      INTEGER J
      DOUBLE PRECISION 
     &   P(NMAX0),
     &   PP(NMAX0),
     &   R(NMAX0),
     &   RR(NMAX0),
     &   Z(NMAX0),
     &   ZZ(NMAX0)
      DOUBLE PRECISION AK,AKDEN,BK,BKDEN,BNRM,BKNUM,
     &   DXNRM,SNRM,XNRM,ZM1NRM,ZNRM
      EXTERNAL SNRM

c ******************************************************************
c - Bi-conjugate gradient method for a sparse system. Adopted from -
c - originally from ``Numerical Recipes''. 
C       USES atimes,asolve,snrm
c history
c 99.04.16 (agli) modified (Princeton Univ.)
c 99.04.22 (agli) modified to make sure that nmax0>=n!!!!!!!!!
c 99.05.17 (agli) modified: now implicit none!!!
c 05.12.14 (BTD) cosmetic cleanup
c end history	
c ******************************************************************
      ITER=0
	
      CALL ATIMES(NMAX,SA,IJA,N,X,R,0)
      DO J=1,N
         R(J)=B(J)-R(J)
         RR(J)=R(J)
      ENDDO

      ZNRM=1.D0
      IF(ITOL.EQ.1)THEN
         BNRM=SNRM(N,B,ITOL)
      ELSEIF(ITOL.EQ.2)THEN
         CALL ASOLVE(NMAX,SA,IJA,N,B,Z,0)
         BNRM=SNRM(N,Z,ITOL)
      ELSEIF(ITOL.EQ.3.OR.ITOL.EQ.4)THEN
         CALL ASOLVE(NMAX,SA,IJA,N,B,Z,0)
         BNRM=SNRM(N,Z,ITOL)
         CALL ASOLVE(NMAX,SA,IJA,N,R,Z,0)
         ZNRM=SNRM(N,Z,ITOL)
      ELSE
         WRITE(0,*)'illegal itol in linbcg'
         STOP
      ENDIF
      CALL ASOLVE(NMAX,SA,IJA,N,R,Z,0)
100   IF(ITER.LE.ITMAX)THEN
         ITER=ITER+1
         ZM1NRM=ZNRM
         CALL ASOLVE(NMAX,SA,IJA,N,RR,ZZ,1)
         BKNUM=0.D0
         DO J=1,N
            BKNUM=BKNUM+Z(J)*RR(J)
         ENDDO
         IF(ITER.EQ.1) THEN
            DO J=1,N
               P(J)=Z(J)
               PP(J)=ZZ(J)
            ENDDO
         ELSE
            BK=BKNUM/BKDEN
            DO J=1,N
               P(J)=BK*P(J)+Z(J)
               PP(J)=BK*PP(J)+ZZ(J)
            ENDDO
         ENDIF
         BKDEN=BKNUM
         CALL ATIMES(NMAX,SA,IJA,N,P,Z,0)
         AKDEN=0.D0
         DO J=1,N
            AKDEN=AKDEN+Z(J)*PP(J)
         ENDDO
         AK=BKNUM/AKDEN
         CALL ATIMES(NMAX,SA,IJA,N,PP,ZZ,1)
         DO J=1,N
            X(J)=X(J)+AK*P(J)
            R(J)=R(J)-AK*Z(J)
            RR(J)=RR(J)-AK*ZZ(J)
         ENDDO
         CALL ASOLVE(NMAX,SA,IJA,N,R,Z,0)
         IF(ITOL.EQ.1.OR.ITOL.EQ.2)THEN
            ZNRM=1.D0
            ERR=SNRM(N,R,ITOL)/BNRM
         ELSEIF(ITOL.EQ.3.OR.ITOL.EQ.4)THEN
            ZNRM=SNRM(N,Z,ITOL)
            IF(ABS(ZM1NRM-ZNRM).GT.EPS*ZNRM)THEN
               DXNRM=ABS(AK)*SNRM(N,P,ITOL)
               ERR=ZNRM/ABS(ZM1NRM-ZNRM)*DXNRM
            ELSE
               ERR=ZNRM/BNRM
               GOTO 100
            ENDIF
            XNRM=SNRM(N,X,ITOL)
            IF(ERR.LE.0.5D0*XNRM) THEN
               ERR=ERR/XNRM
            ELSE
               ERR=ZNRM/BNRM
              GOTO 100
            ENDIF
         ENDIF
c*** diagnostic
c         WRITE(0,*)'bicg ITER=',ITER,' ERR=',ERR
c***
      IF(ERR.GT.TOL)GOTO 100
      ENDIF
      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
c ******************************************************************

      SUBROUTINE ASOLVE(NMAX,SA,IJA,N,B,X,ITRNSP)
      IMPLICIT NONE
      INTEGER NMAX,IJA(NMAX),N,ITRNSP
      DOUBLE PRECISION X(N),B(N),SA(NMAX)
      INTEGER I
      DO I=1,N
         X(I)=B(I)/SA(I)
      ENDDO
      RETURN
      END

c *******************************************************************

      SUBROUTINE ATIMES(NMAX,SA,IJA,N,X,R,ITRNSP)
      IMPLICIT NONE
      INTEGER NMAX,IJA(NMAX),N,ITRNSP
      DOUBLE PRECISION X(N),R(N),SA(NMAX)

C     uses dsprsax,dsprstx

      IF(ITRNSP.EQ.0)THEN
         CALL DSPRSAX(NMAX,SA,IJA,X,R,N)
      ELSE
         CALL DSPRSTX(NMAX,SA,IJA,X,R,N)
      ENDIF
      RETURN
      END
c *******************************************************************
      SUBROUTINE DSPRSAX(NMAX,SA,IJA,X,B,N)
      IMPLICIT NONE
      INTEGER NMAX,IJA(NMAX),N
      DOUBLE PRECISION X(N),B(N),SA(NMAX)

      INTEGER I,K

      IF(IJA(1).NE.N+2)THEN
         WRITE(0,*)'Fatal error: mismatched vector and matrix in sprsax'
         STOP
      ENDIF
      DO I=1,N
         B(I)=SA(I)*X(I)
         DO K=IJA(I),IJA(I+1)-1
            B(I)=B(I)+SA(K)*X(IJA(K))
         ENDDO
      ENDDO
      RETURN
      END
c ********************************************************************
      SUBROUTINE DSPRSTX(NMAX,SA,IJA,X,B,N)
      IMPLICIT NONE
      INTEGER NMAX,IJA(NMAX),N
      DOUBLE PRECISION X(N),B(N),SA(NMAX)

      INTEGER I,J,K

      IF(IJA(1).NE.N+2)THEN
         WRITE(0,*)'fatal error: mismatched vector and matrix in sprstx'
         STOP
      ENDIF
      DO I=1,N
         B(I)=SA(I)*X(I)
      ENDDO
      DO I=1,N
         DO K=IJA(I),IJA(I+1)-1
            J=IJA(K)
            B(J)=B(J)+SA(K)*X(I)
         ENDDO
      ENDDO
      RETURN
      END
c *********************************************************************
 
      FUNCTION SNRM(N,SX,ITOL)
      IMPLICIT NONE
      INTEGER N,ITOL
      DOUBLE PRECISION SNRM,SX(N)	

      INTEGER I,ISAMAX	

      IF(ITOL.LE.3)THEN
         SNRM=0.
         DO I=1,N
            SNRM=SNRM+SX(I)**2
         ENDDO
         SNRM=SQRT(SNRM)
      ELSE
         ISAMAX=1
         DO I=1,N
            IF(ABS(SX(I)).GT.ABS(SX(ISAMAX))) ISAMAX=I
         ENDDO
         SNRM=ABS(SX(ISAMAX))
      ENDIF
      RETURN
      END

c *********************************************************************

      SUBROUTINE DSPRSIN(A,N,THRESH,NMAX,SA,IJA)
      IMPLICIT NONE
      INTEGER N,NMAX,IJA(NMAX)
      DOUBLE PRECISION A(N,N),SA(NMAX),THRESH

      INTEGER I,J,K
	
      DO J=1,N
         SA(J)=A(J,J)
      ENDDO
      IJA(1)=N+2
      K=N+1
      DO I=1,N
         DO J=1,N
            IF(DABS(A(I,J)).GE.THRESH)THEN
               IF(I.NE.J)THEN
                  K=K+1
                  IF(K.GT.NMAX)THEN
                     WRITE(0,*)'fatal error: nmax too small in sprsin'
                     STOP
                  ENDIF
                  SA(K)=A(I,J)
                  IJA(K)=J
               ENDIF
            ENDIF
         ENDDO
         IJA(I+1)=K+1
      ENDDO
      RETURN
      END

