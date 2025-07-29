      SUBROUTINE THERM_COOLING_AFI(I,TI,UI,UIA,UIB,F,UF,UFA,UFB,
     &                             DELTA_U2,NSTATE,NISRF,ISRF_WL,
     &                             ISRF,CABS,AFI_COOLING)

      IMPLICIT NONE

!               therm_cooling_afi_v2

! arguments

      INTEGER I,F,NISRF,NSTATE
      DOUBLE PRECISION 
     &   AFI_COOLING,DELTA_U2,TI,UF,UFA,UFB,UI,UIA,UIB
      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF)

! local variables

      DOUBLE PRECISION 
     &   AFI_COOLING0,AFI_INTRABIN_COOLING,DUF,DUI,
     &   E1,E2,E3,E4,EPS

!-----------------------------------------------------------------------
! subroutine THERM_COOLING_AFI
! given
!     I           = index of initial (upper) bin
!     TI          = temperature (K) associated with initial bin
!     UI          = energy (erg) representative of upper bin
!     UIA,UIB     = lower,upper limits (erg) of upper bin
!     F           = index of final (lower) bin
!     UF          = energy (erg) representative of lower bin
!     UFA,UFB     = lower,upper limits (erg) of lower bin
!     DELTA_U2    = [upper limit of bin 2] - [lower limit of bin 2] (erg)
!     NSTATE      = number of bins used
!     NISRF       = number of wavelengths in tabulated ISRF
!     ISRF_WL[]   = wavelengths (cm) for tabulated ISRF
!     ISRF[]      = c*u_lambda (erg cm-3 s-1) for tabulated ISRF
!     CABS[]      = C_abs (cm2) at wavelengths ISRF_WL

! returns
!     AFI_COOLING = bin-to-bin transition probability (s-1) for
!                   grain in bin I to make transition to bin F
!
! The thermal approximation is used to estimate the transition 
! probability, using either the exact mode distribution or the
! Debye approximation.

! Originally written by Aigen Li, Princeton University
! Modified by B.T. Draine, Princeton University\
! History
! 00.02.18 (AL)  First created
! 00.11.24 (BTD) Cosmetic changes, added comments
! 18.08.14 (BTD) v2 removed A=radius from arg list
! end history
! *******************************************************************

      EPS=1.D-3

! ui > uf: erg; [uia, uib]: boundary of ui;
! [ufa,ufb]: boundary of uf; E1<E2<E3<E4;

      E1=UIA-UFB
      E2=DMIN1((UIA-UFA),(UIB-UFB))
      E3=DMAX1((UIA-UFA),(UIB-UFB))
      E4=UIB-UFA

      DUI=UIB-UIA
      DUF=UFB-UFA

      IF(I.NE.(F+1))THEN

         CALL SIMPSON_THERM_COOLING_RATE(E1,E2,E3,E4,DUI,DUF,TI,
     &                                   NISRF,ISRF_WL,ISRF,CABS,
     &                                   EPS,AFI_COOLING)

         IF(F.NE.1)AFI_COOLING=(AFI_COOLING/(UI-UF))*(UFB-UFA)/(UIB-UIA)
         IF(F.EQ.1)AFI_COOLING=(AFI_COOLING/(UI-UF))*DELTA_U2/(UIB-UIA) 
	
      ELSEIF(I.EQ.(F+1))THEN

         CALL SIMPSON_THERM_COOLING_RATE(E1,E2,E3,E4,DUI,DUF,TI,
     &                                   NISRF,ISRF_WL,ISRF,CABS,
     &                                   EPS,AFI_COOLING0)

         IF(F.NE.1)AFI_COOLING0=(AFI_COOLING0/(UI-UF))*(UFB-UFA)/
     &                          (UIB-UIA)
         IF(F.EQ.1)AFI_COOLING0=(AFI_COOLING0/(UI-UF))*DELTA_U2/
     &                          (UIB-UIA) 

         IF(F.NE.1)THEN
            E1=UIA-UFB
            DUI=UIB-UIA
	    E2=DUI
            CALL SIMPSON_INTRABIN_THERM_COOLING_RATE(E1,E2,DUI,TI,
     &                                               NISRF,ISRF_WL,ISRF,
     &                                               CABS,EPS,
     &                                             AFI_INTRABIN_COOLING)
         ELSEIF(F.EQ.1)THEN	
               AFI_INTRABIN_COOLING=0.D0	
         ENDIF	

         AFI_COOLING=AFI_COOLING0+AFI_INTRABIN_COOLING/(UI-UF)

      ENDIF

      RETURN
      END
