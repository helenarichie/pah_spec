      SUBROUTINE FE_SPEC_HEAT(MODE,NAT,TEMP,EKT,CK)
      IMPLICIT NONE

C Arguments:

      INTEGER MODE
      DOUBLE PRECISION CK,EKT,NAT,TEMP

C Parameters:
C *** Warning: make sure that MXMODES has same value
C     in subroutine sil_temp.f

      INTEGER MXMODES
      PARAMETER(MXMODES=1000000)

C Local variables:

      INTEGER J,NAT1,NAT1OLD,NMODES
      DOUBLE PRECISION FAC,TD3D,X,Y
      DOUBLE PRECISION EMODES(1:MXMODES)
      SAVE EMODES,NMODES
      SAVE NAT1OLD
      DATA NAT1OLD/0/

C***********************************************************************
C
C Given:
C   MODE      = 0 to use mode spectrum computed by subroutine
C                        FE_MODES
C               1 to use 3-D Debye model plus electronic specific heat
C   NAT       = number of atoms in grain
C   TEMP      = T (K)
C
C Returns:
C
C   EKT      = E/kT   , where E(T) = thermal energy at temperature T
C   CK       = C(T)/k , where C(T) = dE/dT
C
C --------------------------------------------------------------
C 
C If MODE=0: We use mode spectrum EMODES(1-NMODES) returned by
C            subroutine SIL_MODES
C We assume that every mode can be approximated as a harmonic
C oscillator:
C
C        N         E_j
C   E = sum [ ------------- ]
C       j=1   exp(E_j/kT)-1
C
C   dE      N     E_j^2 exp(E_j/kT)
C   -- = k sum  -------------------
C   dT     j=1   [exp(E_j/kT)-1]^2
C
C ----------------------------------------------------------------
C
C If MODE=1: We assume specific heat of N atoms to be given by 
C combination 3D Debye model + electronic specific heat.
C
C B.T. Draine, Princeton Univ. Obs., 2011.02.19
C History
C 11.02.19 (BTD) Adapted from sil_spec_heat_v2.f
C End history
C
C***********************************************************************

C Debye temp = 465 K (Desai 1986, JPCRD 15, 967)

      TD3D=465.D0
      FAC=1.D0
      IF(3.*NAT.GT.DBLE(MXMODES))THEN
         FAC=REAL(MXMODES)/(3.D0*NAT)
      ENDIF
      NAT1=INT(FAC*NAT)

      IF(MODE.EQ.0)THEN
         IF(NAT1.NE.NAT1OLD)THEN

C EMODES are mode energies in cm-1

            CALL FE_MODES(MXMODES,NAT1,NMODES,EMODES)
            NAT1OLD=NAT1

c*** NB: this treatment needs to be modified to treat the electronic
c        contribution to the heat capacity

         ENDIF
         EKT=0.D0
         CK=0.D0
         DO J=1,NMODES
            X=1.43877D0*EMODES(J)/TEMP
            IF(X.LT.100.D0)THEN
               Y=EXP(X)
               EKT=EKT+X/(Y-1.D0)
               CK=CK+Y*(X/(Y-1.D0))**2
            ENDIF
         ENDDO
         EKT=EKT/FAC
         CK=CK/FAC
      ELSE

C First compute energy in 3*(NAT-2) vibrational degrees of 
C freedom with 3D Debye model.

         X=TD3D/TEMP
         CALL DEBYE(X,3,EKT,CK)
         EKT=(NAT-2.D0)*3.D0*EKT
         CK=(NAT-2.D0)*3.D0*CK

C Now add electronic specific heat
c gamma   = 4.942*(T/K) mJ/mole*K  (Desai 1986, JPCRD 15, 967)
c gamma/k = 5.944E-4 * (T/K)
c electrons dominate heat content for T < 33.15 K
c vibrations have most of heat for T > 33.15 K

         CK=CK+5.944D-4*NAT*TEMP
         EKT=EKT+2.972D-4*NAT*TEMP

      ENDIF
      RETURN
      END




