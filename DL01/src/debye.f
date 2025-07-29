      SUBROUTINE DEBYE(X,N,EKT,CK)
      IMPLICIT NONE
! Arguments:
      INTEGER N
      DOUBLE PRECISION CK,EKT,X

! Local variables
      DOUBLE PRECISION DEBYEF2,DEBYEF3,X2,X3,XBIG,XSMALL
      EXTERNAL DEBYEF2,DEBYEF3
      DATA XBIG/10./,XSMALL/.01/
!***********************************************************************
! Given:
!     X = T_debye/T
!     N = dimensionality (2 or 3)
! Returns:
!     EKT = energy/kT per degree of freedom
!     CK  = (dE/dT)/k per degree of freedom
!
! B.T. Draine, Princeton University Observatory
! history:
! 99.05.04 BTD first written and tested
! end history
!
!***********************************************************************
      IF(X.GT.XBIG)THEN
         X2=X*X
         IF(N.EQ.2)THEN
            EKT=(2./X2)*(2.*1.20206-X2*EXP(-X)*(1.+2./X+2./X2+
     &                   EXP(-X)*(0.5+0.5/X+0.25/X2)))
         ELSEIF(N.EQ.3)THEN
            X3=X2*X
            EKT=(3./X3)*(6.*1.082323-X3*EXP(-X)*(1.+3./X+6./X2+6./X3+
     &                   EXP(-X)*(0.5+0.75/X+0.75/X2+0.375/X3)))
         ELSE
            STOP'Error in subroutine DEBYE: n is not 2 or 3'
         ENDIF
      ELSEIF(X.LT.XSMALL)THEN
         EKT=1.-0.5*X*N/(N+1.)+N*X*X/(12.*(N+2))
      ELSE
         IF(N.EQ.2)THEN
            CALL DQG32(0.D0,X,DEBYEF2,EKT)
            EKT=2.*EKT/X**2
         ELSEIF(N.EQ.3)THEN
!            write(0,*)'about to call dqg32 with debyef3...'
            CALL DQG32(0.D0,X,DEBYEF3,EKT)
            EKT=3.*EKT/X**3
         ENDIF
      ENDIF

! Have now evaluated energy content E/kT
! Now evaluate specific heat per degree of freedom (1/k)(dE/dT)

      CK=(N+1.)*EKT-N*X/(EXP(X)-1.)
      RETURN      
      END

      DOUBLE PRECISION FUNCTION DEBYEF2(X)
      IMPLICIT NONE
      DOUBLE PRECISION X
      DEBYEF2=X*X/(EXP(X)-1.)
      RETURN
      END

      DOUBLE PRECISION FUNCTION DEBYEF3(X)
      IMPLICIT NONE
      DOUBLE PRECISION X
      DEBYEF3=X**3/(EXP(X)-1.)
      RETURN
      END

