      SUBROUTINE SIL_SPEC_HEAT(MODE,NAT,TEMP,EKT,CK)
      IMPLICIT NONE

! Arguments:

      INTEGER MODE
      DOUBLE PRECISION CK,EKT,NAT,TEMP

! Parameters:
! *** Warning: make sure that MXMODES has same value
!     in subroutine sil_temp.f

      INTEGER MXMODES
      PARAMETER(MXMODES=1000000)

! Local variables:

      INTEGER J,NAT1,NAT1OLD,NMODES
      DOUBLE PRECISION CK2D,CK3D,EKT2D,EKT3D,FAC,F2D,F3D,TD2D,TD3D,X,Y
      DOUBLE PRECISION EMODES(1:MXMODES)
      SAVE EMODES,NMODES
      SAVE NAT1OLD
      DATA NAT1OLD/0/

!***********************************************************************
!
! Given:
!   MODE      = 0 to use mode spectrum computed by subroutine
!                        SIL_MODES
!               1 to use combination of 2-D and 3-D Debye models
!   NAT       = number of atoms in grain
!   TEMP      = T (K)
!
! Returns:
!
!   EKT      = E/kT   , where E(T) = thermal energy at temperature T
!   CK       = C(T)/k , where C(T) = dE/dT
!
! --------------------------------------------------------------
! 
! If MODE=0: We use mode spectrum EMODES(1-NMODES) returned by
!            subroutine SIL_MODES
! We assume that every mode can be approximated as a harmonic
! oscillator:
!
!        N         E_j
!   E = sum [ ------------- ]
!       j=1   exp(E_j/kT)-1
!
!   dE      N     E_j^2 exp(E_j/kT)
!   -- = k sum  -------------------
!   dT     j=1   [exp(E_j/kT)-1]^2
!
! ----------------------------------------------------------------
!
! If MODE=1: We assume specific heat of N atoms to be given by 
! combination of 2D and 3D Debye models.
!
! B.T. Draine, Princeton Univ. Obs., 2000.02.21
! History
! 00.02.21 (BTD) Adapted from pah_spec_heat
! 00.03.31 (BTD) Corrected error -- now check to see whether NAT
!                has changed.
! 08.02.20 (BTD) Modified to allow use with very large grains, with
!                NAT > MXMODES/3
! 22.12.07 (BTD) * Cosmetic changes
! End history
!
!***********************************************************************

! 00.03.27 Fix F2D,F3D=1-F2D,TD2D,TD3D

      F2D=2./3.    ! fraction of modes with 2-D distribution in E
      F3D=1.-F2D   ! fraction of modes with 3-D distribution in E
      TD2D=500.    ! Debye temperature (K) for 2-D mode distribution
      TD3D=1500.   ! Debye temperature (K) for 3-D mode distribution

      FAC=1.D0
      IF(3.*NAT.GT.DBLE(MXMODES))THEN
         FAC=REAL(MXMODES)/(3.D0*NAT)
      ENDIF
      NAT1=INT(FAC*NAT)

      IF(MODE.EQ.0)THEN  ! use mode spectrum computed by SIL_MODES

         IF(NAT1.NE.NAT1OLD)THEN

! EMODES are mode energies in cm-1

            CALL SIL_MODES(MXMODES,NAT1,NMODES,EMODES)
            NAT1OLD=NAT1
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
      ELSE               ! use 2-D and 3-D Debye models for mode spectrum

! First compute energy in 3*F2D*(NAT-2) vibrational degrees of 
! freedom with 2D Debye model.

         X=TD2D/TEMP  ! = T_Debye/T
         CALL DEBYE(X,2,EKT2D,CK2D)
         X=TD3D/TEMP  ! = T_Debye/T
         CALL DEBYE(X,3,EKT3D,CK3D)
         EKT=(NAT-2.D0)*3.D0*(F2D*EKT2D+F3D*EKT3D)
         CK=(NAT-2.D0)*3.D0*(F2D*CK2D+F3D*CK3D)

      ENDIF
      RETURN
      END




