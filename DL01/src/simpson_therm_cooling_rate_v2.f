      SUBROUTINE SIMPSON_THERM_COOLING_RATE(E1,E2,E3,E4,DUI,DUF,T,
     &                                      NISRF,ISRF_WL,ISRF,CABS,
     &                                      EPS,RATE_COOLING)
      IMPLICIT NONE

!               simpson_therm_cooling_rate_v2
! arguments

      INTEGER NISRF
      DOUBLE PRECISION 
     &   DUF,DUI,E1,E2,E3,E4,EPS,RATE_COOLING,T

      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF)

! local variables

      INTEGER K,N
      DOUBLE PRECISION 
     &   ELOW_CUT,EUPP_CUT,H,HC,KB,P,S1,S2,T1,T2,X

! external functions

      DOUBLE PRECISION THERM_COOLING_FUNCT
      EXTERNAL THERM_COOLING_FUNCT

! -----------------------------------------------------------------------
! subroutine SIMPSON_THERM_COOLING_RATE
! given:
!     E1        = Eumin - Elmax (erg)
!     E2        = min(Eumin-Elmin,Eumax-Elmax) (erg)
!     E3        = max(Eumin-Elmin,Eumax-Elmax) (erg)
!     E4        = Eumax - Elmin (erg)
!                where upper bin = [Eumin,Eumax]
!                      lower bin = [Elmin,Elmax]
!     DUI       = Eumax - Eumin (erg)
!     DUF       = Elmax - Elmin (erg)
!     T         = temperature associated with upper bin
!     NISRF     = number of wavelengths for tabulated ISRF
!     ISRF_WL[] = wavelengths (cm) for tabulated ISRF
!     ISRF[]    = c*u_lambda (erg cm-3 s-1) for tabulated ISRF
!     CABS[]    = C_abs (cm2) at wavelengths ISRF_WL
!     EPS       = error tolerance for simpson integration
!
! returns
!     RATE_COOLING = radiated power (erg s-1) for photon energies
!                    in range [E1,E4]
!
! Calculates the radiated power using Simpson's integration method
! using the thermal approximation.
! original code written by Aigen Li, Princeton University
! revised by B.T. Draine, Princeton University
! History
! 00.02.18 (AL)  Original code
! 00.11.24 (BTD) cosmetic changes, added comments
! 18.08.13 (BTD) v2 removed A=radius from argument lists(s)
! end history
! --------------------------------------------------------------- 
! - ---------------------------------------------------------------
! 00.10.28 BTD
!     There does not appear to be any reason to be using Elow_cut and
!     Eupp_cut in calculation of radiative cooling rate, aside from
!     making sure we don't run into trouble with calculation of
!     stimulated emission in function therm_cooling_funct
!     Therefore replace this code redefining E1,E2,E3,E4
!     while simulataneously modifying therm_cooling_funct
!     to avoid trouble with stimulated emission calculation

      HC=6.62607D-27*2.99792D10
      KB=1.38065D-16

      ELOW_CUT=HC/ISRF_WL(NISRF)
      EUPP_CUT=HC/ISRF_WL(1)

      IF(E1.GT.20.*KB*T)THEN
         RATE_COOLING=0.D0
         RETURN
      ENDIF

! --------------------------------------------------------------------  

      N=1
      H=E4-E1

      T1=0.5*H*(THERM_COOLING_FUNCT(E1,T,NISRF,ISRF_WL,ISRF,
     &                              CABS,E1,E2,E3,E4,DUI,DUF)+
     &          THERM_COOLING_FUNCT(E4,T,NISRF,ISRF_WL,ISRF,
     &                              CABS,E1,E2,E3,E4,DUI,DUF))
      S1=T1
10    P=0.0
      DO K=0,N-1
         X=E1+(K+0.5)*H
         P=P+THERM_COOLING_FUNCT(X,T,NISRF,ISRF_WL,ISRF,
     &                           CABS,E1,E2,E3,E4,DUI,DUF)
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

      FUNCTION THERM_COOLING_FUNCT(E,T,NISRF,ISRF_WL,ISRF,
     &                             CABS,E1,E2,E3,E4,DUI,DUF)

      IMPLICIT NONE
      INTEGER NISRF
      DOUBLE PRECISION 
     &   DUI,DUF,E,E1,E2,E3,E4,
     &   T,THERM_COOLING_FUNCT
      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF)

! local variables:

      DOUBLE PRECISION 
     &   C,CRSSCT,GLU_E,H,HC,KB,RADFLD,RADFLD_E,WAVELENGTH

!-----------------------------------------------------------------------
! Given:
!     E = photon energy (erg)
!     T = grain temperature T_u (K) of upper level
!     NISRF = number of wavelengths in tabulated radiation field
!     ISRF_WL(1-NISRF) = wavelengths (cm) for tabulated radiation field
!     ISRF(1-NISRF) = c*u_lambda (erg cm-3 s-1)
!     CABS = C_abs (cm2) for photon energy E
!     E1 = Eumin - Elmax  (erg)
!     E2 = min(Eumin-Elmin,Eumax-Elmax) (erg)
!     E3 = max(Eumin-Elmin,Eumax-Elmax) (erg)
!     E4 = Eumax - Elmin (erg)
!          where upper bin = [Eumin,Eumax]
!                lower bin = [Elmin,Elmax]
!     DUI = Eumax - Eumin (erg)
!     DUF = Elmax - Elmin (erg)
!
! Returns:
!                        DUI*G_lu*E^3*C_abs(E)*[1+(u_E/8*pi)*(h*c/E)^3]
!  THERM_COOLING_FUNCT = ------------------------------------------
!                                 exp(E/kT_u) - 1
!
! where G_lu(E) = function defined by Draine & Li (2001)
!
! Originally written by Aigen Li, Princeton University.
! History:
! 00.11.11 (BTD) Cosmetic changes.  Extensive comments.
! 00.11.24 (BTD) Modified to deal with DUF=0
! 18.08.13 (BTD) v2 removed A = radius from argument list(s)
! end history
!-----------------------------------------------------------------------

      H=6.62607D-27
      C=2.99792D10
      HC=H*C
      KB=1.38065D-16

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

!     evaluate Glu_E
!     dui=width of upper level
!     duf=width of lower level
!     note that Glu_E = G_ul(paper)*dui
!     Note: we assume that DUI > 0 
!     Sanity check:

         IF(DUI.EQ.0.D0)THEN
            WRITE(0,*)'Problem: called THERM_COOLING_FUNCT with DUI=0'
            write(0,*)'E4/hc=',E4/hc
            write(0,*)'E3/hc=',E3/hc
            write(0,*)'E2/hc=',E2/hc
            write(0,*)'E1/hc=',E1/hc
            write(0,*)'E/hc=',E/hc
            write(0,*)'dui,duf=',dui,duf
            STOP
         ENDIF
         IF(E.LT.E1)THEN
            GLU_E=0.D0
         ELSEIF(E.GE.E1 .AND. E.LT.E2)THEN
            GLU_E=(E-E1)/DUF
         ELSEIF(E.GE.E2 .AND. E.LE.E3)THEN
            IF(DUF.GT.0.D0)THEN
               GLU_E=DMIN1(DUI,DUF)/DUF
            ELSE
               GLU_E=1.D0
            ENDIF
         ELSEIF(E.GT.E3 .AND. E.LT.E4)THEN
            GLU_E=(E4-E)/DUF
         ELSE
            GLU_E=0.
         ENDIF

         THERM_COOLING_FUNCT=(GLU_E*(E**3.)/(DEXP(E/(KB*T))-1.D0))*
     &      CRSSCT*(8.*3.1416/(H*HC**2))*
     &      (1.+HC**3/(8.*3.1416*E**3)*RADFLD_E)

      ELSE
         THERM_COOLING_FUNCT=0.D0
      ENDIF
      RETURN
      END
