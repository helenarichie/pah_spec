      SUBROUTINE HEATING_AFI(I,UI,UIA,UIB,F,UF,UFA,UFB,
     &                       NSTATE,NISRF,ISRF_WL,ISRF,
     &                       CABS,AFI_HEATING)

!                      heating_afi_v2

      IMPLICIT NONE

! arguments

      INTEGER I,F,NISRF,NSTATE
      DOUBLE PRECISION AFI_HEATING,UI,UIA,UIB,UF,UFA,UFB
      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF)

! local variables

      DOUBLE PRECISION 
     &   AFI_HEATING0,AFI_HEATING_N,AFI_INTRABIN_HEATING,
     &   DUI,DUF,E1,E2,E3,E4,EC,EPS

!-----------------------------------------------------------------------
! Subroutine HEATING_AFI
! Given:
!     I   = index of initial state
!     UI  = energy (erg) of initial state
!     UIA = lower limit (erg) of initial state
!     UIB = upper limit (erg) of initial state
!     NSTATE = number of energy bins
!     NISRF = number of wavelengths in tabulated ISRF
!     ISRF_WL[] = wavelength (cm) for tabulated ISRF
!     ISRF[] = c*u_lambda (erg cm-3 s-1) for tabulated ISRF
!     CABS[] = C_abs (cm2) at wavelengths of tabulated ISRF
! Returns:
!     AFI_HEATING = contribution to matrix element A_fi (s-1) due
!                   to photon absorptions (heating)
!                 = probability per unit time for transition from
!                   bin I to bin F due to photon absorption

! Originally written by Aigen Li, Princeton University
! History
! 99.09.28 (AL)  original version
! 99.10.21 (AL)  modified
! 00.01.25 (AL)  modified to avoid thermal assumption
! 00.02.07 (AL)  further modified
! 00.02.21 (AL)  include transitions to energies beyond highest bin
! 00.10.26 (AL)  modified to use energy-counting approach
! 00.11.12 (BTD) cosmetic changes; added comments
! 00.11.13 (BTD) rewritten
! 18.08.13 (BTD) v2 remove A from arg list(s)
! end history
!-----------------------------------------------------------------------

      EPS=1.D-3

! E1,E2,E3,E4 are transition energies W1,W2,W3,W4 from Draine & Li (2001)

      E1=UFA-UIB
      E2=DMIN1((UFA-UIA),(UFB-UIB))
      E3=DMAX1((UFA-UIA),(UFB-UIB))
      E4=UFB-UIA
	
      DUI=UIB-UIA
      DUF=UFB-UFA

! calculate heating rate (erg s-1)

      IF(F.NE.NSTATE)THEN
         CALL SIMPSON_HEATING_RATE(E1,E2,E3,E4,DUI,DUF,
     &                             NISRF,ISRF_WL,ISRF,
     &                             CABS,EPS,AFI_HEATING)
      ELSEIF(F.EQ.NSTATE)THEN
         EC=UFA-UIA
         CALL SIMPSON_HEATING_RATE_N(E1,EC,NISRF,ISRF_WL,ISRF,
     &                              CABS,EPS,AFI_HEATING)
      ENDIF
      IF(F.EQ.I+1.AND.I.NE.1)THEN
         CALL SIMPSON_INTRABIN_HEATING_RATE(DUI,NISRF,
     &                                      ISRF_WL,ISRF,CABS,EPS,
     &                                      AFI_INTRABIN_HEATING)
         AFI_HEATING=AFI_HEATING+AFI_INTRABIN_HEATING
      ENDIF

c calculate bin-to-bin transition rate (s-1)

      AFI_HEATING=AFI_HEATING/(UF-UI)

      RETURN
      END
