      SUBROUTINE PAH_SPEC_HEAT(MODE,NC,NH,ND,TEMP,EKT,CK)
      IMPLICIT NONE
c--------------------- pah_spec_heat v2 -------------------------------
C Arguments:

      INTEGER MODE
      DOUBLE PRECISION CK,EKT,NC,ND,NH,TEMP

C Parameters:
C *** Warning: make sure that value of MXMODES agrees with value in
C              subroutine pah_temp.f

      INTEGER MXMODES
      PARAMETER(MXMODES=200000)

C Local variables:

      INTEGER J,NC1,NC1_OLD,ND1,ND1_OLD,NH1,NH1_OLD,NMODES
      DOUBLE PRECISION CK1,EKT1,FAC,NAT,TERM,X,Y
      DOUBLE PRECISION EMODES(1:MXMODES)
      SAVE EMODES,NMODES
      SAVE NC1_OLD,ND1_OLD,NH1_OLD
      DATA NC1_OLD/0/,ND1_OLD/0/,NH1_OLD/0/

C***********************************************************************
C
C Given:
C   MODE      = 0 to use mode spectrum computed by subroutine
C                        PAH_MODES
C               1 to use 2 2-D Debye models for C-C mode spectrum
C                        plus 3 delta functions for C-H modes
C   NC        = number of C atoms (integer)
C   NH        = number of H atoms (integer)
C   ND        = number of D atoms (integer)
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
C            subroutine PAH_MODES
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
C If MODE=1: We assume specific heat of NC-2 atoms to be given by the model of 
C Krumhansl and Brooks (1953): 
C in-plane vibrations characterized by 2-dimensional Debye model
C         with T_Debye = 2500 K
C out-of-plane vibrations characterized by 2-dimensional Debye model
C         T_Debye = 900 K (adjusted downward from 950 K recommended
C                          by Krumhansl & Brooks)
C
C This model closely reproduces measured specific heat of bulk graphite
C for T > 13 K. Experimental data: 
C
C    13 < T < 300 K: DeSorbo & Tyler 1953, J.Chem.Phys. 21, 1660
C    300 < T < 1000 K: data given in Table 1 of Krumhansl & Brooks,
C                      (experimental source unknown)
C
C Assume NH oscillators with hbar*omega/k = 1275 K   (11.3 um)
C        NH     "                           1670 K   ( 8.6 um)
C        NH     "                           4360 K   ( 3.3 um)
C

C B.T. Draine, Princeton Univ. Obs., 2000.02.04
C History
C 00.02.04 (BTD) Created
C 00.02.10 (BTD) Modified to now use list of normal modes returned by
C                pah_modes.f
C 00.02.15 (BTD) Modified to offer option of using either
C                * mode spectrum returned by PAH_MODES
C                * Debye model for C-C modes, plus 3 delta functions
C                  for C-H modes
C 00.03.08 (BTD) Added ND to argument list, and added C-D stretching
C                and bending modes
C 08.02.20 (BTD) Modified to allow treatment of very large numbers of
C                atoms, exceeding MXMODES/3
C 08.02.20 (BTD) v2: modified to have NC,ND,NH as double precision
c                variables
C End history
C
C***********************************************************************

      NAT=NC+ND+NH
      FAC=1.D0
      IF(3.D0*NAT.GT.DBLE(MXMODES))THEN
         FAC=DBLE(MXMODES)/(3.D0*NAT)
      ENDIF
      NC1=NINT(FAC*NC)
      NH1=NINT(FAC*NH)
      ND1=NINT(FAC*ND)
      IF(MODE.EQ.0)THEN
         IF(NC1.NE.NC1_OLD.OR.NH1.NE.NH1_OLD.OR.ND1.NE.ND1_OLD)THEN

C EMODES are mode energies in cm-1

            CALL PAH_MODES(MXMODES,NC1,NH1,ND1,NMODES,EMODES)
            NC1_OLD=NC1
            NH1_OLD=NH1
            ND1_OLD=ND1
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

C First compute energy in NC-2 vibrational degrees of carbon atoms.
C Assume Krumhansl-Brooks model for specific heat of graphite
C Reduce C-C out-of-plane Debye temperature from 950 to 863 K for 
C consistency with model used in PAH_MODES
C 863/1.43877  =  599.8  (cf. 600 cm-1 in PAH_MODES)
C 2500/1.43877 = 1737.6  (cf. 1740 cm-1 in PAH_MODES)
C
C out-of-plane modes:

         X=863.D0/TEMP
         CALL DEBYE(X,2,EKT1,CK1)

C in-plane modes:

         X=2500.D0/TEMP
         CALL DEBYE(X,2,EKT,CK)

C add them together:

         EKT=(NC-2.D0)*(2.D0*EKT+EKT1)
         CK=(NC-2.D0)*(2.D0*CK+CK1)

C Add thermal heat content of H vibrational modes
C 1275/1.43877 =  886.2 cm-1 (cf. 886 cm-1 in PAH_MODES)
C 1670/1.43877 = 1160.7 cm-1 (cf. 1161 cm-1 in PAH_MODES)
C 4360/1.43877 = 3030.4 cm-1 (cf. 3030 cm-1 in PAH_MODES)

         IF(NH.GT.0.D0)THEN
            IF(TEMP.GT.32.D0)THEN
               X=1275.D0/TEMP
               Y=EXP(X)
               TERM=X/(Y-1.)
               EKT1=TERM
               CK1=Y*TERM**2
               IF(TEMP.GT.42.D0)THEN
                  X=1670.D0/TEMP
                  Y=EXP(X)
                  TERM=X/(Y-1.)
                  EKT1=EKT1+TERM
                  CK1=CK1+Y*TERM**2
                  IF(TEMP.GT.109.D0)THEN
                     X=4360./TEMP
                     Y=EXP(X)
                     TERM=X/(Y-1.)
                     EKT1=EKT1+TERM
                     CK1=CK1+Y*TERM**2
                  ENDIF
               ENDIF
               EKT=EKT+NH*EKT1
               CK=CK+NH*CK1
            ENDIF
         ENDIF

C Add thermal heat content of D vibrational modes
C 935.6/1.43877 =  650.3 cm-1 (cf. 886 cm-1 in PAH_MODES)
C 1225.4/1.43877 =  851.7 cm-1 (cf. 1161 cm-1 in PAH_MODES)
C 3200./1.43877 = 2224 cm-1 (cf. 3030 cm-1 for H)

         IF(ND.GT.0.D0)THEN
            IF(TEMP.GT.24.D0)THEN
               X=935.6D0/TEMP
               Y=EXP(X)
               TERM=X/(Y-1.)
               EKT1=TERM
               CK1=Y*TERM**2
               IF(TEMP.GT.31.D0)THEN
                  X=1225.4D0/TEMP
                  Y=EXP(X)
                  TERM=X/(Y-1.)
                  EKT1=EKT1+TERM
                  CK1=CK1+Y*TERM**2
                  IF(TEMP.GT.80.)THEN
                     X=3200./TEMP
                     Y=EXP(X)
                     TERM=X/(Y-1.)
                     EKT1=EKT1+TERM
                     CK1=CK1+Y*TERM**2
                  ENDIF
               ENDIF
               EKT=EKT+ND*EKT1
               CK=CK+NH*CK1
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END

