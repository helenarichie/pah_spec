      SUBROUTINE DBGDL_COOLING_AFI(I,TI,UI,UIA,UIB,F,UF,UFA,UFB,
     &                             NSTATE,NISRF,ISRF_WL,ISRF,CABS,
     &                             AFI_COOLING)

      IMPLICIT NONE

! arguments:

      INTEGER F,I,NISRF,NSTATE
      DOUBLE PRECISION
     &   AFI_COOLING,TI,UI,UIA,UIB,UF,UFA,UFB
      DOUBLE PRECISION
     &   ISRF_WL(NISRF),
     &   ISRF(NISRF),
     &   CABS(NISRF)

      DOUBLE PRECISION EMAX,EMIN,EPS

!-----------------------------------------------------------------------
! subroutine DBDGL_COOLING_AFI
! given
! returns


! *******************************************************************
! - ^^^ Continuous COOLING!!!! Thermal (Debye modes) Method !!! ^^^ -
! - Probability per unit time from state i to f, due to cooling.    -
! - On Feb.18, 2000: thermal (real/Debye modes) treatments.         -
! - On Mar.03, 2000: cooling treated as continous process.          -  
! 18.08.13 (BTD) v2: removed A=radius from arg list
! -------------------------------------------------------------------

      EPS=1.D-3

      IF(I.NE.(F+1))THEN
         AFI_COOLING=0.D0
         RETURN
      ENDIF

      EMIN=0.
      EMAX=UIB

      CALL SIMPSON_DBGDL_COOLING_RATE(EMIN,EMAX,TI,NISRF,ISRF_WL,
     &                                ISRF,CABS,EPS,AFI_COOLING)
      AFI_COOLING=AFI_COOLING/(UI-UF)
	
      RETURN
      END
