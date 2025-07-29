      SUBROUTINE BB_RADIATION_FIELD(TBB,OMEGA,WAVELENGTH,RADFLD)
      IMPLICIT NONE

! arguments

      DOUBLE PRECISION OMEGA,RADFLD,TBB,WAVELENGTH

! local variables

      DOUBLE PRECISION PI

! external functions

      DOUBLE PRECISION PLANCK
      EXTERNAL PLANCK

!-----------------------------------------------------------------------

! BB_RADIATION_FIELD
! given
!     TBB = color temperature (K)
!     OMEGA = dilution factor
!     WAVELENGTH = wavelength (cm)
! returns
!     RADFLD = c*u_lambda (erg cm-3 s-1)
!              for blackbody radiation field cut off at 912 A
! originally written by Aigen Li, Princeton Univ.
! 20.09.10 (BTD) small streamlining
! ---------------------------------------------------------------

      PI=3.14159265D0

      RADFLD=0.
      IF(WAVELENGTH.GE.0.0912D-4)THEN
         RADFLD=OMEGA*4.D0*PI*PLANCK(WAVELENGTH,TBB)
      ENDIF

      RETURN
      END
