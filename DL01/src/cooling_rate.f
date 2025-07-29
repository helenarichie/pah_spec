      SUBROUTINE COOLING_RATE(NISRF,ISRF_WL,
     &                        DLGLAMBDA,CABS,T,RATE_COOLING,SCRW)

      IMPLICIT NONE

! arguments

      INTEGER NISRF
      DOUBLE PRECISION DLGLAMBDA,RATE_COOLING,T
      DOUBLE PRECISION 
     &   CABS(NISRF),
     &   ISRF_WL(NISRF),
     &   SCRW(NISRF)

! local variables

      INTEGER I

! external functions

      DOUBLE PRECISION PLANCK
      EXTERNAL PLANCK

!-----------------------------------------------------------------------
! subroutine COOLING_RATE
! given:
!     NISRF = number of wavelengths in tabulated ISRF
!     ISRF_WL[] = wavelengths (cm) for tabulated ISRF
!                 [assumed to be uniform in lg(lambda)]
!     ISRF[] = c*u_lambda (erg cm-3 s-1) for tabulated ISRF
!     DLGLAMBDA = delta(log10(lambda)) for tabulated ISRF
!     CABS[] = C_abs (cm2) for tabulated ISRF
!     T = grain temperature (K)
! 
! returns:
!     RATE_COOLING = power radiated by grain (erg s-1)
!
! Originally written by Aigen Li, Princeton University
! Modified by B.T. Draine, Princeton University
! History:
! 99.10.11 (AL)  original using Simpson's rule
! 99.10.13 (AL)  replace with trapezoidal integration
! 00.07.19 (AL)  modified
! 00.11.12 (BTD) cosmetic changes, added comments
! 00.11.28 (BTD) add SCRW to arg list, eliminate array RATE_COOLING_WL
! 07.12.11 (BTD) deleted superflous arguments COMPOSITION,A
! end history
!-----------------------------------------------------------------------
!*** diagnostic
!      write(0,*)'cooling_rate ckpt 0, t=',t
!***
      DO I=1,NISRF
	 SCRW(I)=CABS(I)*4.*3.1416*PLANCK(ISRF_WL(I),T)
      ENDDO
!*** diagnostic
!      write(0,*)'cooling_rate ckpt 2, scrw(i)='
!      write(0,fmt='(1pe10.3)')scrw
!      write(0,*)'aobut to call trap_integ with dlglambda=',dlglambda
!***
      CALL TRAP_INTEG(NISRF,ISRF_WL,SCRW,
     &                DLGLAMBDA,RATE_COOLING)

      RETURN
      END
