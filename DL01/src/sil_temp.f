      SUBROUTINE SIL_TEMP(IMODE,ENERGY,NAT,TEMP,EKT,CK)
      IMPLICIT NONE
C
C Arguments:
C
      INTEGER IMODE,NAT
      DOUBLE PRECISION CK,EKT,ENERGY,TEMP
C
C Common:
C
      DOUBLE PRECISION EV1_COM
      COMMON/EVCOM/EV1_COM
      INTEGER MODE_COM,NC_COM,ND_COM,NH_COM,NSIL_COM
      DOUBLE PRECISION EK_COM
      COMMON/GRAINCOM/EK_COM,MODE_COM,NC_COM,ND_COM,NH_COM,NSIL_COM
C
C Local variables:
C
      INTEGER INIT
      DOUBLE PRECISION FZH,FZL,NSIL_DBL,THI,TLO,TOL
      DOUBLE PRECISION T_OLD,TMIN
C
C Arguments of SIL_MODES
C *** Warning: make sure that MXMODES has same value
C     in subroutine sil_spec_heat.f
C
      INTEGER MXMODES
      PARAMETER(MXMODES=1000000)
      INTEGER NMODES
      DOUBLE PRECISION EMODES(1:MXMODES)
C
C External functions:
C
      EXTERNAL FZERO_SIL,ZBRENT_D
      DOUBLE PRECISION FZERO_SIL,ZBRENT_D
C
      DATA INIT/0/,T_OLD/50.D0/
      SAVE EMODES,INIT,T_OLD,TMIN
C
C***********************************************************************
C
C Subroutine SIL_TEMP
C
C Given:
C
C     IMODE   = 0 to find T such that <E(T)>=ENERGY using exact.spec.ht.
C             = 1                                         Debye model
C             = 2 to use T=max[ T s.t. <E(T)>=ENERGY , E_1/(k ln(2)) ]
C                       using exact sp.ht. or Debye model as appropriate
C             = 3 to use T=max[ T s.t. <E(T)>=ENERGY , E_1/(k ln(2)) ]
C                       using Debye model
C             > 3 to use T = E_1/(k ln(2)) for E .le. IMODE*(E_1/10)
C                        T s.t. <E(T)> = ENERGY for E > IMODE*(E_1/10)
C                       using Debye model
C             < 0 to use T = E_1/(k ln(2)) for E .le. EMODE(-IMODE)
C                        T s.t. <E(T)> = ENERGY for E > EMODE(-IMODE)
C
C     ENERGY  = energy (cm-1)
C     NAT     = number of atoms
C
C Returns:
C
C     TEMP    = temperature T (K) of carbon grain
C             = 0 if ENERGY < E_1 = first vibrationally-excited state
C     EK      = <E>/k*T , where <E>=mean thermal energy at temp. T
C     CK      = C/k  , where C(T)=d<E>/dT
C
C Uses subroutine SIL_SPEC_HEAT to relate ENERGY and TEMP
C
C B.T.Draine, Princeton University Observatory
C History:
C 00.03.27 (BTD) adapted from subroutine PAH_TEMP
C 00.03.28 (BTD) further mods
C 00.04.24 (BTD) corrected errors
C 00.04.24 (BTD) modified to bypass call to SIL_MODES when IMODE=1
C 00.04.24 (BTD) modified to set MODE_COM before calling FZERO_SIL
C 00.07.31 (BTD) modified to use more conservative criterion for
C                validity of Debye approx.
C end history
C
C***********************************************************************
C
c*** diagnostic
c      write(0,*)'entered sil_temp with imode=',imode
c***
      EK_COM=1.43877D0*ENERGY
      IF(INIT.EQ.1)THEN
         IF(NSIL_COM.NE.NAT)INIT=0
      ENDIF
      IF(INIT.EQ.0.AND.IMODE.NE.1)THEN
         IF(3*NAT-6.GT.MXMODES)THEN
            WRITE(0,*)'Error in SIL_TEMP: 3*NAT-6 > MXMODES=',MXMODES
            STOP
         ENDIF
         CALL SIL_MODES(MXMODES,NAT,NMODES,EMODES)
         INIT=1
         EV1_COM=EMODES(1)
         TMIN=1.43877*EMODES(1)/LOG(2.D0)
      ENDIF
      NSIL_COM=NAT
      THI=2.*T_OLD
      TLO=0.5*T_OLD
C
C Cannot evaluate temperature for energy less than first excited state
C when using mode spectrum
C
      IF(IMODE.NE.1.AND.ENERGY.LT.0.9999*EMODES(1))THEN
         TEMP=0.
         EKT=0.
         RETURN
      ENDIF
C
C Set MODE_COM=0 to use exact mode distribution for vib. modes
C Set MODE_COM=1 to use 2-D and 3-D Debye models for vib. mode spectrum
C
      IF(IMODE.EQ.0)THEN
         MODE_COM=0
      ELSEIF(IMODE.EQ.1)THEN
         MODE_COM=1
      ELSEIF(IMODE.EQ.2)THEN
         IF(ENERGY.LT.10.*EMODES(1))THEN
            MODE_COM=0
         ELSE
            MODE_COM=1
         ENDIF
      ELSEIF(IMODE.GE.3)THEN
         MODE_COM=1
c      ELSEIF(IMODE.LT.0)THEN
c         IF(ENERGY.LT.10.*EMODES(1))THEN
c            MODE_COM=0
c         ELSE
c            MODE_COM=1
c         ENDIF
      ELSEIF(IMODE.LT.0)THEN
C-------------------------------------------
C modification 00.07.31
C
         IF(ENERGY.LT.NAT*EMODES(1))THEN
            MODE_COM=0
         ELSE
c*** diagnostic
c      write(0,*)'-------------------------- debye'
c***
            MODE_COM=1
         ENDIF
      ELSE
         WRITE(0,*)'Error in SIL_TEMP: IMODE=',IMODE
         STOP
      ENDIF
C
C First carry out crude search to bracket solution
C
 2000 FZL=FZERO_SIL(TLO)
      FZH=FZERO_SIL(THI)
      IF(FZL*FZH.LE.0.D0)THEN
         GOTO 3000
      ELSEIF(FZL.LT.0.D0)THEN
         TLO=4.*TLO
         THI=4.*THI
         GOTO 2000
      ELSE
         TLO=TLO/4.
         THI=THI/4.
         GOTO 2000
      ENDIF
C
C Use ZBRENT to find TEMP
C
 3000 TOL=1.D-3
      TEMP=ZBRENT_D(FZERO_SIL,TLO,THI,TOL)
C
C Save T for use as first guess in next call:
C
      T_OLD=TEMP
      IF(IMODE.EQ.2.OR.
     &   IMODE.EQ.3)THEN
         IF(TEMP.LT.TMIN)TEMP=TMIN
      ELSEIF(IMODE.GE.4)THEN
         IF(ENERGY.LE.IMODE*0.1*EV1_COM)THEN
            TEMP=TMIN
         ENDIF
      ELSEIF(IMODE.LT.0)THEN
         IF(ENERGY.LE.EMODES(-IMODE))THEN
            TEMP=TMIN
         ENDIF
      ENDIF
C
C Compute energy and specific heat for this temperature
C
      NSIL_DBL=NAT
      CALL SIL_SPEC_HEAT(MODE_COM,NSIL_DBL,TEMP,EKT,CK)
      RETURN
      END
      
      DOUBLE PRECISION FUNCTION FZERO_SIL(T)
      IMPLICIT NONE
C
C Arguments:
C
      DOUBLE PRECISION T
C
C Common:
C
      INTEGER MODE_COM,NC_COM,ND_COM,NH_COM,NSIL_COM
      DOUBLE PRECISION EK_COM
      COMMON/GRAINCOM/EK_COM,MODE_COM,NC_COM,ND_COM,NH_COM,NSIL_COM
C
C Local variables:
C
      INTEGER MODE
      DOUBLE PRECISION CK,EKT,NAT_DBL,TPASS

C Given:

C   T = temperature (K)
C
C and, through COMMON/GRAINCOM/ :
C
C   EK_COM  = thermal energy / k
C   MODE    = 0 to use exact mode distribution
C           = 1 to use Debye model for C-C modes
C   NSIL_COM = number of atoms
C
C Returns:
C
C   FZERO_SIL = E(T)/k - EK_COM
C
C           where E(T) = <thermal energy content> of grain at temp T
C           E(T) is computed using a model for the heat capacity
C
      TPASS=T
      MODE=MODE_COM
      NAT_DBL=DBLE(NSIL_COM)
      CALL SIL_SPEC_HEAT(MODE,NAT_DBL,TPASS,EKT,CK)
      FZERO_SIL=EKT*T-EK_COM
      RETURN
      END
