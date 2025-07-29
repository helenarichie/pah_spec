      SUBROUTINE SIMPSON_DBGDL_COOLING_RATE(EMIN,EMAX,T,NISRF,
     &                                      ISRF_WL,ISRF,CABS,EPS,
     &                                      RATE_COOLING)

      IMPLICIT NONE

! arguments:

      INTEGER NISRF
      DOUBLE PRECISION 
     &   EMAX,EMIN,EPS,RATE_COOLING,T
      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF)

! local variables:
      INTEGER K,N
      DOUBLE PRECISION 
     &   ELOW_CUT,EUPP_CUT,H,HC,P,S1,S2,T1,T2,X

! external functions:

      DOUBLE PRECISION DBGDL_COOLING_FUNCT
      EXTERNAL DBGDL_COOLING_FUNCT

!-----------------------------------------------------------------------
! subroutine SIMPSON_DBGDL_COOLING_RATE
! given
!     EMIN
!     EMAX
!     T
!     NISRF
!     ISRF_WL[]
!     ISRF[]
!     CABS[]
!     EPS
! returns:
!     RATE_COOLING =
!
! Originally written by Aigen Li, Princeton University
! History
! 00.02.18 (AL)  modified
! 00.11.26 (BTD) cosmetic changes, added comments
! 18.08.13 (BTD) v2 removed A=radius from arg list
! end history
!-----------------------------------------------------------------------

      HC=6.62607D-27*2.99792D10
      ELOW_CUT=HC/ISRF_WL(NISRF)
      EUPP_CUT=HC/ISRF_WL(1)

      IF(EMIN.LT.ELOW_CUT)EMIN=ELOW_CUT
      IF(EMIN.GT.EUPP_CUT)EMIN=EUPP_CUT

      IF(EMAX.LT.ELOW_CUT)EMAX=ELOW_CUT
      IF(EMAX.GT.EUPP_CUT)EMAX=EUPP_CUT
      
      IF(EMIN.EQ.EMAX)THEN
         RATE_COOLING=0.D0
         RETURN
      ENDIF

      N=1
      H=EMAX-EMIN
      T1=0.5*H*(DBGDL_COOLING_FUNCT(EMIN,T,NISRF,ISRF_WL,ISRF,CABS)+
     &          DBGDL_COOLING_FUNCT(EMAX,T,NISRF,ISRF_WL,ISRF,CABS))
      S1=T1
 10   P=0.0
      DO K=0,N-1
         X=EMIN+(K+0.5)*H
         P=P+DBGDL_COOLING_FUNCT(X,T,NISRF,ISRF_WL,ISRF,CABS)
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

      FUNCTION DBGDL_COOLING_FUNCT(E,T,NISRF,ISRF_WL,ISRF,CABS)

      IMPLICIT NONE
      INTEGER NISRF
      DOUBLE PRECISION
     &   DBGDL_COOLING_FUNCT,E,T
      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF)

! local variables

      DOUBLE PRECISION 
     &   C,CRSSCT,H,HC,KB,PI,RADFLD,RADFLD_E,WAVELENGTH

!-----------------------------------------------------------------------
! function DBGDL_COOLING_FUNCT
! given
!     E
!     A
!     T
!     NISRF
!     ISRF_WL
!     ISRF
!     CABS
! returns
!     DBGDL_COOLING_FUNCT =

! Originally written by Aigen Li, Princeton University
! History
! 00.02.18 (AL)  modified
! 00.03.02 (AL)  modified
! 00.11.26 (BTD) cosmetic changed, added comments
! 18.08.13 (BTD) v2 removed A=radius from arg list
! end history
!-----------------------------------------------------------------------

      H=6.62607D-27
      C=2.99792D10
      HC=H*C
      KB=1.38065D-16
      PI=3.14159D0

      IF(E.NE.0.)THEN

         WAVELENGTH=H*C/E
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

         RADFLD_E=RADFLD*WAVELENGTH**2/(H*C*C)

         DBGDL_COOLING_FUNCT=(E**3.)/(DEXP(E/(KB*T))-1.D0)*
     &                       CRSSCT*(8.*PI/(H*HC**2))*
     &                       (1.+HC**3/(8.*PI*E**3)*RADFLD_E)

      ELSEIF(E.EQ.0.)THEN

         DBGDL_COOLING_FUNCT=0.D0

      ENDIF
      RETURN
      END
