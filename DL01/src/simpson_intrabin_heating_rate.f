      SUBROUTINE SIMPSON_INTRABIN_HEATING_RATE(DUI,NISRF,
     &                                         ISRF_WL,ISRF,CABS,EPS,
     &                                         RATE_HEATING)

      IMPLICIT NONE

c arguments

      INTEGER NISRF
      DOUBLE PRECISION DUI,EPS,RATE_HEATING
      DOUBLE PRECISION CABS(NISRF),ISRF(NISRF),ISRF_WL(NISRF)

c local variables

      INTEGER K,N
      DOUBLE PRECISION 
     &   E1,E2,ELOW_CUT,EUPP_CUT,H,HC,P,S1,S2,T1,T2,X

c external functions

      DOUBLE PRECISION INTRABIN_HEATING_FUNCT
      EXTERNAL INTRABIN_HEATING_FUNCT
c-----------------------------------------------------------------------
c subroutine SIMPSON_INTRABIN_HEATING_RATE
c given:
c     DUI = width of lower bin (erg)
c     NISRF = number of wavelengths for tabulated ISRF
c     ISRF_WL[] = wavelengths (cm) for tabulated ISRF
c     ISRF[] = c*u_lambda (erg cm-3 s-1) for tabulated ISRF
c     CABS[] = C_abs (cm2) at wavelengths ISRF_WL
c     EPS = fractional error tolerance

c returns
c                     DUI       E
c     RATE_HEATING = \int (1 - ---)*C_abs(E)*c*u_E*dE
c                      0       DUI
c
c                DUI       E                              dlambda
c             = \int (1 - ---)*C_abs(E)*c*u_lambda*lambda*-------
c                 0       DUI                              lambda
c
c since E*u_E = lambda*u_lambda 
c         u_E = (lambda/E)*u_lambda
c         u_E = (lambd
c Originally written by Aigen Li, Princeton University
c Rewritten by B.T. Draine, Princeton University
c History
c 00.11.07 (AL)  First written
c 00.11.13 (BTD) Rewritten.
c end history
c --------------------------------------------------------------- 

      HC=6.62607D-27*2.99792D10

      IF(DUI.LE.0.D0)THEN
         RATE_HEATING=0.D0
         RETURN
      ENDIF

      ELOW_CUT=HC/ISRF_WL(NISRF)
      EUPP_CUT=HC/ISRF_WL(1)

      IF(ELOW_CUT.GT.DUI)THEN
         RATE_HEATING=0.D0
         RETURN
      ELSE
         E1=ELOW_CUT
         E2=DUI
         IF(DUI.GT.EUPP_CUT)E2=EUPP_CUT
      ENDIF

c --------------------------------------------------------------------  
      N=1
      H=E2-E1
      T1=0.5*H*(
     &    INTRABIN_HEATING_FUNCT(E1,DUI,NISRF,ISRF_WL,ISRF,CABS)+
     &    INTRABIN_HEATING_FUNCT(E2,DUI,NISRF,ISRF_WL,ISRF,CABS))
      S1=T1
 10   P=0.0
      DO K=0,N-1
         X=E1+(K+0.5)*H
         P=P+INTRABIN_HEATING_FUNCT(X,DUI,NISRF,ISRF_WL,ISRF,CABS)
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
      RATE_HEATING=S2

      RETURN
      END

c***********************************************************************

      FUNCTION INTRABIN_HEATING_FUNCT(E,DUI,NISRF,ISRF_WL,ISRF,CABS)

      IMPLICIT NONE

c arguments

      DOUBLE PRECISION INTRABIN_HEATING_FUNCT

      INTEGER NISRF
      DOUBLE PRECISION DUI,E
      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF)

c local variables

      DOUBLE PRECISION C,CRADFLD_E,CRSSCT,H,RADFLD,WAVELENGTH

c-----------------------------------------------------------------------
c function INTRABIN_HEATING_FUNCT
c given
c     E = photon energy (erg)
c     DUI = width of lower bin (erg)
c     NISRF = number of wavelengths in tabulated ISRF
c     ISRF_WL[] = wavelengths (cm) for tabulated ISRF
c     ISRF[] = c*u_lambda (erg cm-3 s-1) for tabulated ISRF
c     CABS[] = C_abs (cm2) at wavelengths ISRF_WL[]
c returns
c                                   E
c     INTRABIN_HEATING_FUNCT = (1- ---)*C_abs*c*u_E
c                                  DUI
c where
c               lambda**2
c       c*u_E = --------- * c*u_lambda
c                  h*c
c
c Originally written by Aigen Li, Princeton University
c Modified by B.T. Draine, Princeton University
c history:
c 00.11.13 (BTD) Modified, added comments
c end history
c-----------------------------------------------------------------------

      H=6.62607D-27
      C=2.99792D10
	
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

         CRADFLD_E=RADFLD*WAVELENGTH**2/(H*C)
         INTRABIN_HEATING_FUNCT=(1.-E/DUI)*CRSSCT*CRADFLD_E

      ELSEIF(E.EQ.0.)THEN

         INTRABIN_HEATING_FUNCT=0.D0

      ENDIF
	
      RETURN
      END
