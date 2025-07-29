      SUBROUTINE SIMPSON_INTRABIN_THERM_COOLING_RATE(E1,E2,DUI,TI,
     &                                               NISRF,ISRF_WL,ISRF,
     &                                               CABS,EPS,
     &                                               RATE_COOLING)

!             simpson_intrabin_therm_cooling_rate_v2

      IMPLICIT NONE

! arguments:

      INTEGER NISRF
      DOUBLE PRECISION 
     &   DUI,E1,E2,EPS,RATE_COOLING,TI
      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF)
      INTEGER K,N
      DOUBLE PRECISION 
     &   ELOW_CUT,EUPP_CUT,H,HC,KB,P,S1,S2,T1,T2,X

! external functions:

      DOUBLE PRECISION INTRABIN_THERM_COOLING_FUNCT
      EXTERNAL INTRABIN_THERM_COOLING_FUNCT

!-----------------------------------------------------------------------
! subroutine SIMPSON_INTRABIN_THERM_COOLING_RATE
! given:
!     E1           =
!     E2           =
!     A            =
!     DUI          =
!     TI           =
!     NISRF        =
!     ISRF_WL[]    =
!     ISRF[]       =
!     CABS[]       =
!     EPS          =

! returns
!     RATE_COOLING = intrabin cooling rate (erg s-1) calculated using
!                    Simpson's rule integration

! Originally written by Aigen Li, Princeton University
! History:
! 00.02.18 (AL)  modified
! 00.10.28 (BTD) modified
! 00.11.26 (BTD) cosmetic changes, added comments
! 18.08.13 (BTD) v2 removed A=radius from arg list(s)
! end history
!-----------------------------------------------------------------------

      HC=6.62607D-27*2.99792D10
      KB=1.38065D-16

      ELOW_CUT=HC/ISRF_WL(NISRF)
      EUPP_CUT=HC/ISRF_WL(1)

      IF(E1.GT.20.*KB*TI)THEN
         RATE_COOLING=0.D0
         RETURN
      ENDIF
C --------------------------------------------------------------------  

      N=1
      H=E2-E1

      T1=0.5*H*(INTRABIN_THERM_COOLING_FUNCT(E1,DUI,TI,NISRF,
     &                                       ISRF_WL,ISRF,CABS)+
     &          INTRABIN_THERM_COOLING_FUNCT(E2,DUI,TI,NISRF,
     &                                       ISRF_WL,ISRF,CABS))

      S1=T1
 10   P=0.0
      DO K=0,N-1
         X=E1+(K+0.5)*H
         P=P+INTRABIN_THERM_COOLING_FUNCT(X,DUI,TI,NISRF,
     &                                    ISRF_WL,ISRF,CABS)
      ENDDO
      T2=(T1+H*P)/2.0
      S2=(4.0*T2-T1)/3.0
      IF((DABS((S2-S1)/S1)).GE.EPS)THEN
         T1=T2
         N=N+N
         H=H/2.0
         S1=S2
         GOTO 10
      ENDIF
      RATE_COOLING=S2

      RETURN
      END

!***********************************************************************

      FUNCTION INTRABIN_THERM_COOLING_FUNCT(E,DUI,TI,NISRF,
     &                                      ISRF_WL,ISRF,CABS)

      IMPLICIT NONE

! arguments:

      INTEGER NISRF
      DOUBLE PRECISION 
     &   DUI,E,INTRABIN_THERM_COOLING_FUNCT,TI
      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF)

! local variables:

      DOUBLE PRECISION 
     &   C,CRSSCT,H,HC,KB,PI,RADFLD,RADFLD_E,WAVELENGTH

!-----------------------------------------------------------------------
! function INTRABIN_THERM_COOLING_FUNCT
! given
!     E
!     A
!     DUI
!     TI
!     NISRF
!     ISRF_WL[]
!     ISRF[]
!     CABS[]

! returns
!     INTRABIN_THERM_COOLING_FUNCT =

! Originally written by Aigen Li, Princeton University
! History
! 00.02.18 (AL)  modified
! 00.02.28 (AL)  modified
! 00.10.28 (BTD) modified
! 00.11.26 (BTD) cosmetic changes, added comments
! 18.08.13 (BTD) v2 removed A=radius from arg list
! end history
!-----------------------------------------------------------------------

      H=6.62607D-27
      C=2.99792D10
      HC=H*C
      KB=1.38065D-16
      PI=3.14159D0

      IF(E.NE.0.)THEN

         WAVELENGTH=HC/E

         IF(WAVELENGTH.GE.ISRF_WL(1) .AND. 
     &      WAVELENGTH.LE.ISRF_WL(NISRF))THEN
            CALL PARAB_INTERPL(ISRF_WL,CABS,NISRF,WAVELENGTH,CRSSCT)
            CALL PARAB_INTERPL(ISRF_WL,ISRF,NISRF,WAVELENGTH,RADFLD)
         ELSEIF(WAVELENGTH.LT.ISRF_WL(1))THEN
            CRSSCT=CABS(1)
            RADFLD=0.
         ELSEIF(WAVELENGTH.GT.ISRF_WL(NISRF))THEN
            CRSSCT=CABS(NISRF)*(ISRF_WL(NISRF)/WAVELENGTH)**2
            RADFLD=0.
         ENDIF	

         RADFLD_E=RADFLD*WAVELENGTH**2/(HC*C)	

         INTRABIN_THERM_COOLING_FUNCT=(E**3)*CRSSCT*
     &                                (8.*PI/(H*HC**2))*
     &                                (1.-E/DUI)*
     &                                (1.+HC**3/(8.*PI*E**3)*RADFLD_E)/
     &                                (EXP(E/(KB*TI))-1.)

      ELSEIF(E.EQ.0.)THEN

         INTRABIN_THERM_COOLING_FUNCT=0.D0

      ENDIF

      RETURN
      END
