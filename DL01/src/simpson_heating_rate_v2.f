      SUBROUTINE SIMPSON_HEATING_RATE(E1,E2,E3,E4,DUI,DUF,
     &                                NISRF,ISRF_WL,ISRF,
     &                                CABS,EPS,RATE_HEATING)

      IMPLICIT NONE

! arguments:

      INTEGER NISRF
      DOUBLE PRECISION 
     &   DUI,DUF,E1,E2,E3,E4,EPS,RATE_HEATING
      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF)

! local variables

        INTEGER N,K
        DOUBLE PRECISION ELOW_CUT,EUPP_CUT,H,HC,P,S1,S2,T1,T2,X

! external functions

        DOUBLE PRECISION HEATING_FUNCT
	EXTERNAL HEATING_FUNCT

!-----------------------------------------------------------------------
! Subroutine SIMPSON_HEATING_RATE
! given:
!     E1 = U_fmin - U_imax (erg)               = W1 in Draine & Li (2001)
!     E2 = U_fmin - U_imin (erg) (if DUF < DUI) = W2 in DL2001
!     E3 = U_fmax - U_imax (erg) (if DUF < DUI) = W3 in DL2001
!     E4 = U_fmax - U_imin (erg)                = W4 in DL2001
!     DUI = width of initial (lower) bin (erg)
!     DUF = width of final (upper) bin (erg)
!     A  = radius (cm)
!     NISRF = number of wavelengths in tabulated ISRF
!     ISRF_WL[] = wavelengths (cm) for tabulated ISRF
!               (assumed uniform in lg(lambda))
!     ISRF[] = c*u_lambda (erg cm-3 s-1) for tabulated ISRF
!     CABS[] = C_abs (cm2) at wavelengths of tabulated ISRF
!     EPS  = error tolerance for Simpson integration

! returns

!     RATE_HEATING = rate of energy absorbed in transitions from
!                    bin I (lower) to bin F (upper)

! Originally written by Aigen Li, Princeton University
! History:
! 99.11.06 (AL)  Modified
! 00.01.25 (AL)  Modified
! 00.02.07 (AL)  Modified
! 00.11.12 (BTD) cosmetic changes; added comments
!                some streamlining
! 18.08.13 (BTD) v2 remove A from arg list
! end history
!-----------------------------------------------------------------------
      HC=6.62607D-27*2.99793D10

      IF(E1.EQ.E4)THEN
         RATE_HEATING=0.D0
         RETURN
      ENDIF	

      ELOW_CUT=HC/ISRF_WL(NISRF)
      EUPP_CUT=HC/ISRF_WL(1)

      IF(E1.LT.ELOW_CUT)E1=ELOW_CUT
      IF(E1.GT.EUPP_CUT)E1=EUPP_CUT

      IF(E2.LT.ELOW_CUT)E2=ELOW_CUT
      IF(E2.GT.EUPP_CUT)E2=EUPP_CUT

      IF(E3.LT.ELOW_CUT)E3=ELOW_CUT
      IF(E3.GT.EUPP_CUT)E3=EUPP_CUT

      IF(E4.LT.ELOW_CUT)E4=ELOW_CUT
      IF(E4.GT.EUPP_CUT)E4=EUPP_CUT
	
! --------------------------------------------------------------------  
      N=1
      H=E4-E1
      T1=0.5*H*(HEATING_FUNCT(E1,NISRF,ISRF_WL,ISRF,CABS,
     &                        E1,E2,E3,E4,DUI,DUF)+
     &          HEATING_FUNCT(E4,NISRF,ISRF_WL,ISRF,CABS,
     &                        E1,E2,E3,E4,DUI,DUF))
      S1=T1
 10   P=0.0
      DO K=0,N-1
         X=E1+(K+0.5)*H
         P=P+HEATING_FUNCT(X,NISRF,ISRF_WL,ISRF,CABS,
     &                     E1,E2,E3,E4,DUI,DUF)
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

!***********************************************************************

      FUNCTION HEATING_FUNCT(E,NISRF,ISRF_WL,ISRF,CABS,
     &                       E1,E2,E3,E4,DUI,DUF)

      IMPLICIT NONE
      DOUBLE PRECISION HEATING_FUNCT

! arguments

      INTEGER NISRF
      DOUBLE PRECISION DUI,DUF,E,E1,E2,E3,E4
      DOUBLE PRECISION 
     &   ISRF_WL(NISRF),
     &   ISRF(NISRF),
     &   CABS(NISRF)

! local variables

      DOUBLE PRECISION 
     &   C,CRSSCT,GUL_E,H,RADFLD,RADFLD_E,WAVELENGTH

!-----------------------------------------------------------------------
! function HEATING_FUNCT
! given:
!     E = photon energy (erg)
!     NISRF = number of wavelengths in tabulated ISRF
!     ISRF_WL[] = wavelengths (cm) for tabulated ISRF
!                 (assumed uniform in lg(lambda))
!     ISRF[] = c*u_lambda (erg cm-3 s-1) for tabulated ISRF
!     CABS[] = C_abs (cm2) at wavelengths ISRF_WL[]
!     E1,E2,E3,E4 = transition energies W1,W2,W3,W4 (erg)
!                   defined by Draine & Li (2001)
!     DUI = width (erg) of initial (lower) bin
!     DUF = width (erg) of final (upper) bin
!
! returns:
!                     c*dUf
!     HEATING_FUNCT = -----*G_ul(E)*C_abs(E)*u_E
!                     Uf-Ui        
!
!                      dUf                  lambda**2
!                   = -----*G_ul(E)*Cabs(E)*---------*c*u_lambda
!                     Uf-Ui                   h*c**2
!
!  so that T_ul = \int HEATING_FUNCT dE
!
! Originally written by Aigen Li, Princeton University
! Modified by B.T. Draine, Princeton University
! History
! 00.11.12 (BTD) cosmetic changes, added comments, streamline
! end history
!-----------------------------------------------------------------------

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

         RADFLD_E=RADFLD*WAVELENGTH**2/(H*C)

         IF(E1.NE.E2)THEN

            IF(E.LT.E1)THEN
               GUL_E=0.D0
            ELSEIF(E.GE.E1 .AND. E.LT.E2)THEN
               GUL_E=(E-E1)/DUI
            ELSEIF(E.GE.E2 .AND. E.LE.E3)THEN
               IF(E2.EQ.E3)THEN
                  GUL_E=1.
               ELSEIF(E2.NE.E3)THEN
                  GUL_E=DMIN1(DUI,DUF)/DUI
               ENDIF
            ELSEIF(E.GT.E3 .AND. E.LT.E4)THEN
               GUL_E=(E4-E)/DUI
            ELSEIF(E.GE.E4)THEN
               GUL_E=0.D0
            ENDIF

            HEATING_FUNCT=GUL_E*CRSSCT*RADFLD_E

         ELSEIF(E1.EQ.E2)THEN

            IF(E.LT.E1 .OR. E.GE.E4)THEN
               HEATING_FUNCT=0.D0
            ELSEIF(E.GE.E1 .AND. E.LT.E4)THEN
               HEATING_FUNCT=CRSSCT*RADFLD_E
            ENDIF
         ENDIF
		
      ELSEIF(E.EQ.0.)THEN

         HEATING_FUNCT=0.D0

      ENDIF

      RETURN
      END

      SUBROUTINE SIMPSON_HEATING_RATE_N(E1,EC,NISRF,ISRF_WL,ISRF,
     &                                  CABS,EPS,RATE_HEATING)

      IMPLICIT NONE

! arguments

      INTEGER NISRF
      DOUBLE PRECISION E1,EC,EPS,RATE_HEATING
      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF)

! local variables

      INTEGER K,N
      DOUBLE PRECISION 
     &   ELOW_CUT,EMIN,EUPP_CUT,H,HC,P,S1,S2,T1,T2,X

! external functions

      DOUBLE PRECISION HEATING_FUNCT_N
      EXTERNAL HEATING_FUNCT_N

!-----------------------------------------------------------------------
! subroutine SIMPSON_HEATING_RATE_N
! given
!     E1 = UNA - UIB (erg)
!     EC = UNA - UIA (erg)
!          where UNA = lower bound of bin N
!                UIA = lower bound of bin I
!                UIB = upper bound of bin I
!     NISRF = number of wavelengths in tabulated ISRF
!     ISRF_WL[] = wavelengths (cm) for tabulated ISRF
!     ISRF[] = c*u_lambda (erg cm-3 s-1) for tabulated ISRF
!     CABS[] = C_abs (cm2) at wavelengths ISRF_WL[]
!     EPS = fractional error tolerance
! returns
!     RATE_HEATING = energy absorption rate (erg s-1) assigned to
!                    transition from initial energy bin to highest 
!                    energy bin, including real
!                    transitions to levels beyond highest energy bin
!
!
! Original version written by Aigen Li, Princeton University
! Rewritten by B.T. Draine, Princeton University
! History
! 00.11.13 (BTD) rewrite
! end history
!-----------------------------------------------------------------------
      HC=6.62607D-27*2.99792D10

! Elow_cut: the lower limit of IR photon (cut at 8000 um) (00.02.22).	
! Eupp_cut: the upper limit of UV photon (cut at 13.6eV) (00.02.21).

      ELOW_CUT=HC/ISRF_WL(NISRF)
      EUPP_CUT=HC/ISRF_WL(1)

      EMIN=E1
      IF(EMIN.LT.ELOW_CUT)EMIN=ELOW_CUT
	 
      IF(EMIN.GE.EUPP_CUT)THEN
         RATE_HEATING=0.D0
         RETURN
      ENDIF

      N=1
      H=EUPP_CUT-EMIN
      T1=0.5*H*(HEATING_FUNCT_N(EMIN,E1,EC,NISRF,ISRF_WL,ISRF,CABS)+
     &          HEATING_FUNCT_N(EUPP_CUT,E1,EC,NISRF,ISRF_WL,ISRF,CABS))
      S1=T1
 10   P=0.0
      DO K=0,N-1
         X=EMIN+(K+0.5)*H
         P=P+HEATING_FUNCT_N(X,E1,EC,NISRF,ISRF_WL,ISRF,CABS)
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

!***********************************************************************

      FUNCTION HEATING_FUNCT_N(E,E1,EC,NISRF,ISRF_WL,ISRF,CABS)

      IMPLICIT NONE

! arguments

      DOUBLE PRECISION HEATING_FUNCT_N

      INTEGER NISRF
      DOUBLE PRECISION E,E1,EC
      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF)

! local variables

      DOUBLE PRECISION 
     &   C,CRADFLD_E,CRSSCT,H,RADFLD,WAVELENGTH

!-----------------------------------------------------------------------
! function HEATING_FUNCT_N
! given
!     E = photon energy (erg)
!     E1 = UNA - UIB (erg)
!     EC = UNA - UIA (erg)
!          where UNA = lower bound of bin N
!                UIA = lower bound of bin I
!                UIB = lower bound of bin I
!     NISRF = number of wavelengths in tabulated ISRF
!     ISRF_WL[] = wavelengths (cm) for tabulated ISRF
!     ISRF[] = c*u_lambda (erg cm-3 s-1) for tabulated ISRF
!     CABS[] = C_abs (cm2) at wavelengths ISRF_WL[]
! returns
!     HEATING_FUNCT_N = 
!        
!-----------------------------------------------------------------------
! Originally written by Aigen Li, Princeton University
! Modified by B.T. Draine, Princeton University
! history
! 00.11.13 (BTD) minor modifications, cosmetic changes, added comments
! end history
!-----------------------------------------------------------------------

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

         IF(E.LT.EC)THEN
            HEATING_FUNCT_N=CRSSCT*CRADFLD_E*(E-E1)/(EC-E1)
         ELSE
            HEATING_FUNCT_N=CRSSCT*CRADFLD_E
         ENDIF

      ELSEIF(E.EQ.0.)THEN
         HEATING_FUNCT_N=0.D0	
      ENDIF
	
      RETURN
      END
