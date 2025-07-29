      SUBROUTINE PAH_TEMP(IMODE,ENERGY,NC,NH,ND,TEMP,EKT,CK)
      IMPLICIT NONE
c---------------- pah_temp v2 ------------------------------------
C Arguments:

      INTEGER IMODE,NC,ND,NH
      DOUBLE PRECISION CK,EKT,ENERGY,TEMP

C Common:

      DOUBLE PRECISION EV1_COM
      COMMON/EVCOM/EV1_COM
      INTEGER MODE_COM,NC_COM,ND_COM,NH_COM,NSIL_COM
      DOUBLE PRECISION EK_COM
      COMMON/GRAINCOM/EK_COM,MODE_COM,NC_COM,NH_COM,ND_COM,NSIL_COM

C Local variables:

      INTEGER INIT,NC_INT,ND_INT,NH_INT
      DOUBLE PRECISION FZH,FZL,NC_DBL,ND_DBL,NH_DBL,THI,TLO,TOL
      DOUBLE PRECISION T_OLD,TMIN

C Arguments of PAH_MODES
C *** Warning: make sure that value of MXMODES agrees with value in
C              subroutine pah_spec_heat.f

      INTEGER MXMODES
      PARAMETER(MXMODES=1000000)
      INTEGER NMODES
      DOUBLE PRECISION EMODES(1:MXMODES)

C External functions:

      EXTERNAL FZERO_PAH,ZBRENT_D
      DOUBLE PRECISION FZERO_PAH,ZBRENT_D

      DATA INIT/0/,T_OLD/50.D0/
      SAVE EMODES,INIT,T_OLD,TMIN

C***********************************************************************
C
C Subroutine PAH_TEMP
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
C     NC      = number of C atoms
C     NH      = number of H atoms
C     ND      = number of D 
C
C Returns:
C
C     TEMP    = temperature T (K) of carbon grain
C             = 0 if ENERGY < E_1 = first vibrationally-excited state
C     EK      = <E>/k*T , where <E>=mean thermal energy at temp. T
C     CK      = C/k  , where C(T)=d<E>/dT
C
C Uses subroutine PAH_SPEC_HEAT to relate ENERGY and TEMP
C
C B.T.Draine, Princeton University Observatory
C History:
C 99.05.10 (BTD) Created this module to invert E(T) provided by
C                subroutine DEBYE
C 99.05.26 (BTD) modified to include vibrational modes of NH H atoms
C 99.05.26 (BTD) modified to use separate subroutine PAH_SPEC_HEAT
C                to characterize thermal properties of grain
C                Removed COMMON/DEBYECOM/...
C                Added COMMON/GRAINCOM/NC_COM,NH_COM
C 00.02.10 (BTD) modified to use new version of pah_spec_heat,
C                which uses mode spectrum produced by pah_modes
C 00.03.08 (BTD) modified to add ND as argument, and as part of
C                argument list for PAH_SPEC_HEAT
C 00.03.27 (BTD) added IMODE to argument list of PAH_TEMP
C                and now allow selection of method for estimating
C                temperature
C                added MODE to COMMON/GRAINCOM/
C 00.03.27 (BTD) further mods to experiment with using
C                T = E_1/(k ln 2)  for E .le. IMODE*E_1
C                T such that <E>=E for E .gt. IMODE*E_1
C 00.03.28 (BTD) further modifications
C 00.04.24 (BTD) modified to bypass call to PAH_MODES when IMODE=1
C 00.04.24 (BTD) modified to set MODE_COM before calling FZERO_PAH
C 00.07.19 (BTD) modified to require E > NC*EMODES(1) as condition
C                to use Debye heat capacity calculation
c 00.11.21 (BTD) modified criterion for returning T=0 for low
c                energies
c 06.05.23 (BTD) changed TOL=1.D-3 -> TOL=1.D-5
C end history
C
C***********************************************************************
C
      EK_COM=1.43877D0*ENERGY
      IF(INIT.EQ.1)THEN
         IF(NC_COM.NE.NC)INIT=0
         IF(NH_COM.NE.NH)INIT=0
         IF(ND_COM.NE.ND)INIT=0
      ENDIF
C
C If IMODE=1, do not need to call PAH_MODES to initialize EMODES
C             but otherwise need to do so
C
      IF(INIT.EQ.0.AND.IMODE.NE.1)THEN
         IF(3*(NC+NH+ND)-6.GT.MXMODES)THEN
            WRITE(0,*)'Error in PAH_TEMP: 3*(NC+NH+ND)-6 > MXMODES=',
     &                MXMODES
            STOP
         ENDIF
         CALL PAH_MODES(MXMODES,NC,NH,ND,NMODES,EMODES)
         INIT=1
         EV1_COM=EMODES(1)
         TMIN=1.43877*EMODES(1)/LOG(2.D0)
      ENDIF
      NC_COM=NC
      NH_COM=NH
      ND_COM=ND
      THI=2.*T_OLD
      TLO=0.5*T_OLD
C
C Cannot evaluate temperature for energy less than first excited state
C when using mode spectrum
C
c--------------------------------------------------------
c modified 00.11.21
c      IF(INIT.NE.1.AND.ENERGY.LT.0.9999*EMODES(1))THEN
      IF(ENERGY.LE.0.0D0)THEN
c--------------------------------------------------------
         TEMP=0.
         EKT=0.
         RETURN
      ENDIF
C
C Set MODE_COM=0 to use exact mode distribution for C-C modes
C Set MODE_COM=1 to use 2-D Debye models for C-C mode spectrum
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
      ELSEIF(IMODE.LT.0)THEN
C-------------------------------------------
C modification 00.07.19
C
C         IF(ENERGY.LT.10.*EMODES(1))THEN
         IF(ENERGY.LT.NC*EMODES(1))THEN
C-------------------------------------------
            MODE_COM=0
         ELSE
c*** diagnostic
c      write(0,*)'------------------------ debye'
c***
            MODE_COM=1
         ENDIF
      ELSE
         WRITE(0,*)'Error in PAH_TEMP: IMODE=',IMODE
         STOP
      ENDIF
C
C First carry out crude search to bracket solution
C
 2000 FZL=FZERO_PAH(TLO)
      FZH=FZERO_PAH(THI)
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

c-----------------------------------------------------------------------
c 06.05.23 (BTD) change TOL=1.D-3 -> TOL=1.D-5
c 3000 TOL=1.D-3
 3000 TOL=1.D-5
c-----------------------------------------------------------------------
      TEMP=ZBRENT_D(FZERO_PAH,TLO,THI,TOL)
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

      NC_DBL=DBLE(NC)
      ND_DBL=DBLE(ND)
      NH_DBL=DBLE(NH)
      CALL PAH_SPEC_HEAT(MODE_COM,NC_DBL,NH_DBL,ND_DBL,TEMP,EKT,CK)
      RETURN
      END
      
      DOUBLE PRECISION FUNCTION FZERO_PAH(T)
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
      COMMON/GRAINCOM/EK_COM,MODE_COM,NC_COM,NH_COM,ND_COM,NSIL_COM
C
C Local variables:
C
      INTEGER MODE
      DOUBLE PRECISION CK,EKT,NC_DBL,ND_DBL,NH_DBL,TPASS
C
C Given:
C
C   T = temperature (K)
C
C and, through COMMON/GRAINCOM/ :
C
C   EK_COM  = thermal energy / k
C   MODE_COM= 0 to use exact mode distribution
C           = 1 to use Debye model for C-C modes
C   NC_COM  = number of C atoms with vibrational degrees of freedom
C   NH_COM  = number of H atoms
C   ND_COM  = number of D atoms
C
C Returns:
C
C   FZERO_PAH = E(T)/k - EK_COM
C
C           where E(T) = <thermal energy content> of grain at temp T
C           E(T) is computed using a model for the heat capacity
C
      TPASS=T
      MODE=MODE_COM
      NC_DBL=NC_COM
      NH_DBL=NH_COM
      ND_DBL=ND_COM
      CALL PAH_SPEC_HEAT(MODE,NC_DBL,NH_DBL,ND_DBL,TPASS,EKT,CK)
      FZERO_PAH=EKT*T-EK_COM
      RETURN
      END













