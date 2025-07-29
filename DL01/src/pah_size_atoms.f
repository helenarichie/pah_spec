      SUBROUTINE PAH_SIZE_ATOMS(A,NC,NH,HTOC)
      IMPLICIT NONE
      DOUBLE PRECISION A,HTOC,NC,NH

c We assume radius to be related to NC by
c NC = 460*(a/1.e-7cm)**2   <--> carbon density 2.20 g cm-3
C
c Noting that coronene = C_24H_12 is pericondensed with H/C=0.5
c                        C_52H_18 is pericondensed with H/C=0.346
c                        C_96H_24 is pericondensed with H/C=0.25
C
c We will assume H/C = 0.5 for N_C < 25
c                H/C = 0.5/(N_C/25)^{1/2} for 25 < N_C < 100
c                H/C = 0.25 for 100 < N_C 
C
c history
c 07.12.11 (BTD) changed assumed density of carbon from 
c                2.24 to 2.20 g cm-3
C 08.02.20 (BTD) changed NC,NH from integer to DOUBLE PRECISION
c                to avoid risk of overflow for a > 1.3e-5 cm
c 20.06.16 (BTD) changed assume density of carbon from
c                2.20 gm cm-3 to 2.00 g cm-3
c end history
c ---------------------------------------------------------------
! old code
!      NC=460.D21*A**3
! 20.06.16 (BTD) change to rho_C = 2.0 g cm-3
      NC=418.D21*A**3

      IF(NC .LE. 25.D0)THEN
         NH=0.5*NC
      ELSEIF(NC.GT.25.D0.AND.NC.LE.100.D0)THEN
         NH=0.5*SQRT(25.*NC) 
      ELSEIF(NC.GT.100.D0)THEN
         NH=0.25*NC
      ENDIF

c*** btd change 00.11.20

      HTOC=NH/NC

c***
      RETURN
      END
