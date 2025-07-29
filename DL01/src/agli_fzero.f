      FUNCTION AGLI_FZERO(NISRF,ISRF_WL,DLGLAMBDA,CABS,
     &                    RATE_HEATING,T,SCRW)
      IMPLICIT NONE

c arguments:

      DOUBLE PRECISION AGLI_FZERO
      INTEGER NISRF
      DOUBLE PRECISION DLGLAMBDA,RATE_HEATING,T
      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF_WL(NISRF),
     &   SCRW(NISRF)

c local variables

      DOUBLE PRECISION RATE_COOLING

c-----------------------------------------------------------------------
c function AGLI_FZERO
c given:
c     COMPOSITION = 'PAH' or 'sil'
c     A = radius (cm)
c     NISRF = number of wavelengths in tabulated ISRF
c     ISRF_WL[] = wavelength (cm) for tabulate ISRF
c                 (assumed to be uniform in lg(lambda))
c     ISRF[] = c*u_lambda (erg cm-3 s-1) for tabulated ISRF
c     DLGLAMBDA = delta(log10(lambda)) for tabulated ISRF
c     CABS[] = C_abs (cm2) for wavelengths ISRF_WL
c     RATE_HEATING = heating rate (erg s-1) for grain in ISRF
c returns
c     AGLI_FZERO = (Pcool-Pheat)/(Pcool+Pheat)
c                  where Pcool = cooling rate for grain
c                        Pheat = heating rate for grain
c
c Originally written by Aigen Li, Princeton University
c History
c 00.11.12 (BTD) cosmetic changes, comments added
c 00.11.28 (BTD) add SCRW to argument list, and to arg list of
c                COOLING_RATE
c 07.12.11 (BTD) remove superfluous COMPOSITION,A from argument list
c                of COOLING RATE
c                removed superfluous COMPOSITION,A from argument list
c                of AGLI_FZERO
c end history
c-----------------------------------------------------------------------
      CALL COOLING_RATE(NISRF,ISRF_WL,DLGLAMBDA,CABS,
     &                  T,RATE_COOLING,SCRW)

      AGLI_FZERO=(RATE_COOLING-RATE_HEATING)/
     &           (RATE_COOLING+RATE_HEATING)

      RETURN
      END
