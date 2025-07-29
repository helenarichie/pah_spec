      SUBROUTINE MMP_RADIATION_FIELD(WAVELENGTH,RADFLD)
      IMPLICIT NONE

! arguments

      DOUBLE PRECISION WAVELENGTH,RADFLD

! local variables

      DOUBLE PRECISION 
     & C,H,K,W2,W3,W4,PI,
     & TBB2,TBB3,TBB4,WAVELENGTHUM

      DOUBLE PRECISION PLANCK
      EXTERNAL PLANCK

c-----------------------------------------------------------------------
! subroutine MMP_RADIATION
! given
!     WAVELENGTH = wavelength (cm)
! returns
!     RADFLD = c*u_lambda (erg cm-3 s-1)
!              for the interstellar radiation field in the solar
!              neighborhood
!              (Mathis et al. 1983; Mezger et al. 1982)

! Originally written by Aigen Li, Princeton University
! Modified by B.T. Draine, Princeton University
! history
! 00.11.13 (BTD) cosmetic changes, added comments
! 08.02.02 (BTD) corrected numerical value of OMEGA3:
!                was OMEGA3=1.D-13 but should be 1.65D-13
!                changed TBB5 from 2.9 to 2.725
! 20.09.10 (BTD) changed OMEGA -> W
!                changed W4=4.D-13 -> 7.D-13
!                remove CMB
! end history
!-----------------------------------------------------------------------
      H=6.62607D-27
      C=2.99792D10
      K=1.38065D-16
      PI=4.D0*ATAN(1.D0)

      W2=1.D-14
      TBB2=7500.

      W3=1.65D-13
      TBB3=4000.

      W4=7.D-13
      TBB4=3000.

      WAVELENGTHUM=WAVELENGTH*1.D4

      IF(WAVELENGTHUM .LT. 0.0912)THEN
         RADFLD=0.
      ELSEIF(WAVELENGTHUM.GE.0.0912 .AND. WAVELENGTHUM.LT.0.11)THEN
         RADFLD=38.57*WAVELENGTHUM**3.4172*1.D4
      ELSEIF(WAVELENGTHUM.GE.0.11 .AND. WAVELENGTHUM.LT.0.134)THEN
         RADFLD=0.02045*1.D4
      ELSEIF(WAVELENGTHUM.GE.0.134 .AND. WAVELENGTHUM.LE.0.246)THEN
         RADFLD=7.115D-4*WAVELENGTHUM**(-1.6678)*1.D4
      ELSEIF(WAVELENGTHUM.GT.0.246)THEN
         RADFLD=4.*PI*(W2*PLANCK(WAVELENGTH,TBB2)+
     &                 W3*PLANCK(WAVELENGTH,TBB3)+
     &                 W4*PLANCK(WAVELENGTH,TBB4))
      ENDIF

      RETURN
      END
