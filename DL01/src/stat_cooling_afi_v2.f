      SUBROUTINE STAT_COOLING_AFI(I,UI,UIA,UIB,LNGUI,F,UF,UFA,UFB,
     &                            LNGUF,NSTATE,NISRF,ISRF_WL,ISRF,
     &                            CABS,AFI_COOLING)

!                stat_cooling_afi_v2

      IMPLICIT NONE

! arguments

      INTEGER I,F,NSTATE,NISRF
      DOUBLE PRECISION 
     &   AFI_COOLING,LNGUF,LNGUI,
     &   UF,UFA,UFB,UI,UIA,UIB
      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF)

! local variables

      DOUBLE PRECISION 
     &   AFI_COOLING0,AFI_INTRABIN_COOLING,
     &   C,DUI,DUF,E1,E2,E3,E4,EPS,H

!-----------------------------------------------------------------------
! subroutine STAT_COOLING_AFI
! given:
!     I,F         = index of initial, final state
!     UI,UF       = midbin energy (erg) of initial, final state
!     UIA,UIB     = lower, upper energy (erg) of bin i
!     UFA,UFB     = lower, upper energy (erg) of bin f
!     LNGUI,LNGUF = ln(g_i), ln(g_f) where g_i,g_f = degeneracy of
!                   bin i,f
!     NSTATE      = number of energy bins
!     NISRF       = number of wavelength in tabulated ISRF
!     ISRF_WL[]   = wavelength (cm) for tabulated ISRF
!     ISRF[]      = c*u_lambda (erg cm-3 s-1) for tabulated ISRF
!     CABS[]      = C_abs (cm2) at wavelengths isrf_wl[]

! returns:
!     AFI_COOLING = transition rate (s-1) from bin i to bin f
!                   due to radiative cooling
!                   calculated using "exact-statistical" method
!
! this routine originally written by Aigen Li, Princeton Univ.
! further modified by B.T. Draine, Princeton University
! History:
! 99.09.28 (AL)  Original version
! 00.02.18 (AL)  Statistical version created
! 00.11.20 (BTD) cosmetic changes, added comments
! 18.08.13 (BTD) v2 removed A = radius from arg list
! end history
!-----------------------------------------------------------------------
      H=6.62607D-27
      C=2.99792D10
      EPS=1.D-3

! E1<E2<E3<E4;

      E1=UIA-UFB
      E2=DMIN1((UIA-UFA),(UIB-UFB))
      E3=DMAX1((UIA-UFA),(UIB-UFB))
      E4=UIB-UFA
	
      DUI=UIB-UIA
      DUF=UFB-UFA

      IF(I.NE.(F+1))THEN

         CALL SIMPSON_STAT_COOLING_RATE(E1,E2,E3,E4,DUI,DUF,NISRF,
     &                                  ISRF_WL,ISRF,CABS,EPS,
     &                                  AFI_COOLING)
         IF((LNGUF-LNGUI).LT.-1.D300 .OR. (LNGUF-LNGUI).GT.1.D300)THEN
            AFI_COOLING=0.D0
            WRITE(0,FMT='(A,2I6)')'stat_cooling_afi_v2 ckpt 1: I,F=',I,F

! 00.10.29 btd fix:

         ELSEIF((LNGUF-LNGUI).GE.-1.D300 .AND. 
     &          (LNGUF-LNGUI).LE.1.D300)THEN
            AFI_COOLING=AFI_COOLING/(UI-UF)*DEXP(LNGUF-LNGUI)
         ENDIF

      ELSEIF(I.EQ.(F+1))THEN

         CALL SIMPSON_STAT_COOLING_RATE(E1,E2,E3,E4,DUI,DUF,NISRF,
     &                                  ISRF_WL,ISRF,CABS,EPS,
     &                                  AFI_COOLING0)
	 IF((LNGUF-LNGUI).LT.-1.D300 .OR. (LNGUF-LNGUI).GT.1.D300)THEN
	    AFI_COOLING0=0.D0
            WRITE(0,FMT='(A,2I6)')'stat_cooling_afi_v2 ckpt 2: I,F=',I,F
         ELSEIF((LNGUF-LNGUI).GE.-1.D300 .AND. 
     &          (LNGUF-LNGUI).LE.1.D300)THEN
            AFI_COOLING0=AFI_COOLING0/(UI-UF)*DEXP(LNGUF-LNGUI)	
         ENDIF

         E1=UIA-UFB
	 DUI=UIB-UIA
         E2=DUI
	 CALL SIMPSON_INTRABIN_STAT_COOLING_RATE(E1,E2,DUI,NISRF,
     &                                           ISRF_WL,ISRF,CABS,EPS,
     &                                           AFI_INTRABIN_COOLING)

         AFI_COOLING=AFI_COOLING0+AFI_INTRABIN_COOLING/(UI-UF)

      ENDIF

      RETURN
      END
