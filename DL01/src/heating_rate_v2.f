      SUBROUTINE HEATING_RATE(NISRF,ISRF_WL,ISRF,DLGLAMBDA,
     &                        CABS,RATE_HEATING,RATE_PHOTONS,SCRW)

!    ------------------ heating_rate_v2 -----------------
      IMPLICIT NONE

! arguments

      INTEGER NISRF
      DOUBLE PRECISION DLGLAMBDA,RATE_HEATING,RATE_PHOTONS
      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF),
     &   SCRW(NISRF)

! local variables

      INTEGER I
      DOUBLE PRECISION HC

!-----------------------------------------------------------------------
! subroutine HEATING_RATE
! given:

!     NISRF = number of wavelengths in tabulated ISRF
!     ISRF_WL[] = wavelengths (cm), uniform in lg(lambda)
!     ISRF[] = c*u_lambda (erg cm-3 s-1) at wavelengths ISRF_WL
!     DLGLAMBDA = delta(log_10(lambda)) for wavelengths ISRF_WL
!     CABS[] = C_abs (cm2) at wavelengths ISRF_WL[]
!     SCRW[] = scratch array

! returns:

!     RATE_HEATING = heating rate (erg s-1)

! Originally written by Aigen Li, Princeton University
! Modified by B.T. Draine, Princeton University
! History: 
! 99.10.06 (AL)  Originally written
! 99.10.13 (AL)  Simpson's rule integration replaced by trapezoidal
! 00.07.19 (AL)  Modified
! 00.11.12 (BTD) Modified; delete MODE option from argument list
!                (not needed since we only calculate heating rates)
!                cosmetic changes, added comments
!                [Why are we using trapezoidal integration instead of
!                simpson's rule???]
! 07.12.11 (BTD) eliminated superfluous arguments:
!                A,COMPOSITION
! 16.04.26 (BTD) v2: add calculation of photon absorption rate
! end history
!-----------------------------------------------------------------------
!*** diagnostic
!      write(0,*)'heating_rate_v2 ckpt 1: NISRF=',NISRF
!      write(0,fmt='(a,1pe10.3,a,e10.3)')'     isrf_wl(300)=',
!     &   isrf_wl(300),' isrf(300)=',isrf(300)
!      write(0,*)'  isrf_wl     isrf     cabs'
!***
      DO I=1,NISRF
         SCRW(I)=CABS(I)*ISRF(I)
!*** diagnostic
!         write(0,fmt='(1p3e10.3)')isrf_wl(i),isrf(i),cabs(i)
!***
      ENDDO

! calculate heating rate

      CALL TRAP_INTEG(NISRF,ISRF_WL,SCRW,
     &                DLGLAMBDA,RATE_HEATING)

! calculate photoabsorption rate

      DO I=1,NISRF
         SCRW(I)=CABS(I)*ISRF(I)*ISRF_WL(I)
      ENDDO
      CALL TRAP_INTEG(NISRF,ISRF_WL,SCRW,
     &                DLGLAMBDA,RATE_PHOTONS)
      HC=1.98645D-16
      RATE_PHOTONS=RATE_PHOTONS/HC

!*** diagnostic
      write(0,fmt='(a,1pe10.3,a,1pe10.3,a,a,1pe10.3,a)')
     &   'heating_rate_v2 ckpt 2: rate_heating=',rate_heating,
     &   ' erg/s; photon rate=',rate_photons,' s-1',
     &   ' <hnu>=',(rate_heating/(1.602e-12*rate_photons)),' eV'
!***
      RETURN
      END
