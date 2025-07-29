      PROGRAM VSG_STAT_THERM_V13
      IMPLICIT NONE
!------------------- vsg_stat_therm_v13 -----------------------------
      INTEGER NCASEMX
      PARAMETER(NCASEMX=42)
      CHARACTER
     &   CASEDESCRIPT*25,
     &   CDESCR(NCASEMX)*79,
     &   COMPOSITION*25,
     &   GEOM*4,
     &   METHOD*5,
     &   RADIATION*7,
     &   RADDESCRIPT*30

! parameters

      INTEGER NDIME,NDIMMODES,NISRF0,NMAX0,NSIZE0,NSTATEMX
      PARAMETER(NDIME=1000000)

! generous allocation of memory for NDIMMODES -- up to 333,000 atoms

      PARAMETER(NDIMMODES=1000000)
!      PARAMETER(NSTATEMX=1002)
      PARAMETER(NSTATEMX=501)
      PARAMETER(NISRF0=2500)

! 121004 (BTD) allow for doubling of number of sizes
! 170606 (BTD) back to original
! 200627 (BTD) back to 200
!      PARAMETER(NSIZE0=100)

      PARAMETER(NSIZE0=200)
      
! Derived parameters

      PARAMETER(NMAX0=NSTATEMX*NSTATEMX)

! Local variables

      INTEGER 
     &   I,ICAS,ICASE,ICGREMT2,IMODE,ISHAP,ISIZE,IWAV,IWRITE,
     &   J,J1,J2,JEFF,JJ,JNEW,JS,JSZ,JZONE,
     &   NISRF,NMAX,NMID,NSTATE,NSIZE,NSTATE1,NSTATEZ1,NZONE

      INTEGER
     &   IJA(NMAX0)

      REAL
     &   TAUNH,TAUVNH,WAVEUM

      DOUBLE PRECISION 
     &   AGE,AMIN,AMAX,AV,AVMAX,BOVERA,C,CHI,
     &   DA,DAV,DEDT,DELTA_LGA,DHRDAV,DLGLAMBDA,DTAU,DTOH,
     &   EBAR,EEQ,FAC,H,HC,HR,ICMB,
     &   LAMBDAMIN,LAMBDAMAX,LOCATION,LUMINOSITY,
     &   OMEGA,PI,PORO,PSUM,RADIUS_STAR,RATE_HEATING,
     &   SUM,SUMD,TCMB,TEFF,TSTAR,UCMB,UISRF,UMMP,USTARLIGHT,
     &   WGT,WAVE_EFF,ZBC,ZBPASS,DELTAMIN,DELTAMAX,LAMBD_ULAMBD

      DOUBLE PRECISION
     &   AWAV(1:NISRF0),
     &   A(1:NSIZE0),
     &   ABSORPTION(1:NISRF0),
     &   AMATRIX(1:NSTATEMX,1:NSTATEMX),
     &   AMATRIX1(1:NSTATEMX,1:NSTATEMX),
     &   BGD(1:NSTATEMX,1:NSTATEMX),
     &   CABS(1:NISRF0),
     &   COOLING_TIME(NSTATEMX),
     &   CSCA(1:NISRF0),
     &   DELTA_T(1:NSTATEMX),
     &   E(0:NDIME),
     &   EMISSION(1:NISRF0),
     &   EMISSIONZ(1:NISRF0),   ! new
     &   EMODES(1:NDIMMODES),
     &   ISRF(1:NISRF0),
     &   ISRF0(1:NISRF0),
     &   ISRFZ(1:NISRF0),
     &   ISRF_WL(1:NISRF0),
     &   LNG(0:NDIME),
     &   LNGSUM(0:NDIME),
     &   LNGU(1:NSTATEMX),
     &   PCUMULT(1:NSTATEMX),
     &   PSTATE(1:NSTATEMX),
     &   PSTATEZ(1:NSTATEMX),
     &   PSTATE1(1:NSTATEMX),
     &   RAD_TIME(1:NSTATEMX),
     &   RHS(1:NSTATEMX),
     &   SA(1:NMAX0),
     &   SCRW(1:NISRF0),
     &   T(1:NSTATEMX),
     &   TZ(1:NSTATEMX),
     &   U(1:NSTATEMX),
     &   UA(1:NSTATEMX),
     &   UB(1:NSTATEMX),
     &   UZ(1:NSTATEMX),
     &   UAZ(1:NSTATEMX),
     &   UBZ(1:NSTATEMX)

! External functions

      DOUBLE PRECISION PLANCK
      EXTERNAL PLANCK
!-----------------------------------------------------------------------
! Program VSG_STAT_THERM
! Purpose: to compute the energy distribution function for grains which
! are small enough that quantized heating broadens the distribution.
! Program can calculate energy distributions for following compositions:
!
!  * Carbonaceous grains ("PAH"-like for small grains, graphite-like
!                         for large grains)
!  * Graphite grains
!  * Silicate grains
!  * Fe grains
! [* amorphous carbon grains --- not yet...]
!
! For carbonaceous grains, optical cross sections depend on charge state
! so it is therefore necessary to specify the charge state.
! Method is as described in Draine & Li (2001), with choice of
! following methods:

!    'dbdis' : "thermal-discrete approximation"

!    'dbcon' : "thermal-continuous approximation"

!    'stati' : "exact-statistical"

!    'dbgdl' : another version of Debye-thermal-continuous, using slightly
!              different but faster-to-evaluate expression for the downward
!              transition rates

! For carbonaceous grains, we use optical
! cross sections as described in Li & Draine (2001). 

! History: originally developed by Aigen Li, Princeton University.
! Further modified by B.T. Draine, Princeton University.
! History:
! 99.09.28 (AL)  Enthalpy range treated self-consistently
! 99.10.18 (AL)  modified enthalpy bins, 
!                treat both silicates and graphite
! 99.10.29 (AL)  adopt Joe's size distribution
! 99.11.27 (AL)  adopt new charge calculation
! 99.11.30 (AL)  modifed to put PAH charge distribution calculation
!                in subroutine and valid for CNM+WNM+WIM with MMP ISRF
! 00.01.25 (AL)  substantially modified -- no thermal approximation
!                use Einstein coefficients
! 00.02.18 (AL)  consider both statistical and thermal treatments in
!                single code.
!                For thermal treatment, two cases:
!                1. derive T from real vibrational modes
!                2. derive T from Debye model for modes
! 00.03.03 (AL)  Add continuous cooling method
! 00.03.08 (AL)  Add option of deuterated PAHs
! 00.03.20 (AL)  Separate Debye-continuous cooling case into two
!                methods: 
!                a. dbcon A[u-1,u]=sum_A[l,u]*(Eu-El)/(Eu-Eu-1)
!                b. dbgdl A[u-1,u]=
!                   int_0^Eu{4.*pi*B(Tu)*Cabs*dlambda}/(Eu-Eu-1)
! 00.07.19 (AL)  Just calculate P(T) and emission for one size
! 00.10.18 (AL)  Modify to allow treatment of large sizes (a > 200A).
! 00.11.11 (BTD) Cosmetic changes to make code more readable
!                with extensive comments added.
!                Modify to input all optional parameters from input
!                file 'vsg_stat_therm.par'
! 00.11.15 (BTD) continuing editing
! 00.11.24 (BTD) modify for separate input files:
!                vsg_isrf.par to specify radiation field
!                vsg_meth.par to specify computational method
!                vsg_size.par to specify size and composition
! 00.11.28 (BTD) major modification to improve handling of memory
!                allocation by passing reusable memory arrays to
!                subroutines instead of having independent definitions
!                In particular, this allowed elimination of arrays
!                E_ORI(0:NDIME),LNG_ORI(0:NDIME),LNGSUM_ORI(0:NDIME).
!                For NDIME=10^6, this reduces memory use by 24 MBytes!
!                Also modify in order to have all adjustable
!                dimensioning done in main program and passed to
!                subroutines and functions.
! 05.10.30 (BTD) added OPEN(UNIT=66,FILE='vsg_stat_therm.pout')
!                added OPEN(UNIT=77,FILE='vsg_stat_therm.iout')
!                [to support WRITE(66... and WRITE(77... statements
!                in subroutine VSG_TD_EMISSION
!                added CLOSE(66) and CLOSE(77) statements
! 05.10.31 (BTD) added OPEN(UNIT=70,FILE='vsg_stat_therm.dpdtlib.out')
!                added CLOSE(70)
! 06.06.05 (BTD) added 'PAD' as new composition option:
!                PAD = PAH with D/(H+D)=0.25
! 07.12.10 (BTD) changed type of ZG to INTEGER
! 08.03.08 (BTD) moved initialization of NSTATE1 so that it is initialized
!                to NSTATE-1 before each call of VSG_TD_EMISSION
!                (this is necessary because when called for PAHs, PAH 
!                neutrals can enter steady-state domain (for which 
!                VSG_TD_EMISSION resets NSTATE1 to NSTATE1=1)
!                before this applies to PAH ions.  Therefore need to reset
!                NSTATE1 to NSTATE-1 before calling VSG_TD_EMISSION for
!                ions.
! 11.02.19 (BTD) v4
!                * Adapted from v3
!                * Modified to allow treatment of Fe grains
! 12.10.04 (BTD) v5
!                * Adapted from v4
!                * Modified to allow radiation field = M31BLGE
!                  (spectrum of M31 bulge)
!                * Modified to allow doubling of number of grain sizes
!                  can now use Delta(log10(a))=0.025
! 12.10.11 (BTD) * Correctly normalize bulge radiation field heating
!                  rate 
!                  U=1 now gives heating rate for MW3.1_60 dust mix
!                  that is same as MMP83
! 12.11.11 (BTD) * add calculation of starlight energy density,
!                  with output to standard out
! 14.02.17 (BTD) v6
!                * added ICASE,ISHAP as variables
!                * add line to read ISHAP from vsg_size.par
!                * add ICASE,ISHAP to argument list of VSG_TD_EMISSION
!                * change COMPOSITION*3  -> COMPOSITION*8
!                * eliminate variable COMPOSITION0
!                * eliminate variable CHARGE_STATUS
!                  not needed now that we have three compositions
!                  gra_PAHi, gra_PAHn, and gra_PAHs
!                * eliminate ZG
!                * remove ZG from arg list of
!                  * vsg_td_emission_v6
!                  * transition_matrix_v3
! 16.02.08 (BTD) v7 modified to handle new cases
!                  ICASE=19,20,21,22
!                * modified to do ICASE and COMPOSITION sanity check
! 16.02.20 (BTD) * modified to handle ICASE=23-33
! 17.12.25 (BTD) v8 updated NCASEMX to NCASEMX=42
! 18.02.28 (BTD) v9 modified to use new ICASE numbers
! 18.06.08 (BTD) * added new cases 209-218, 294-298
! 18.08.11 (BTD) v10 modified to use new cases DH18sil_Px.xx_x.xx_Fex.xx
! 18.10.08 (BTD) * added new cases for f_8=0.7, P=0,0.1,0.2,0.3
! 19.01.21 (BTD) v11 ensure support for ICASE=569 (D16 graphite, MG2 EMT)
! 19.01.31 (BTD) * change to new scheme f8=0.4,0.5,..., 1.1,1.2
! 19.08.01 (BTD) v12 change from DH18 astrosil to DH19 Astrodust
!                    DH19Ad_Px.xx_x_xx
!                Poro from 0 to 0.9 in steps 0.1
!                f_Fe from 0 to 0.5 in steps 0.1
! 19.08.07 (BTD) * add support for BRUZUAL_CHARLOT_RF
! 20.03.06 (BTD) * modify comments to update to DH20Ad_Px.xx_x.xx models
! 20.05.10 (BTD) * remove some unused arrays
!                * introduce new array EMISSIONZ
!                * eliminate PDR option (do this with SLAB)
!                * for SLAB option, use EMISSIONZ to augment
!                  EMISSION = average emission per grain, averaged over
!                  total thickness AVMAX
! 20.05.20 (BTD) * modify to write file 'isrf.out' with adopted radiation
!                  field
! 20.06.09 (BTD) * modify to write out file extcurv.out (using UNIT=13)
!                  in cases where extinction matters ('SLAB' or 'LIRG')
! 20.06.29 (BTD) * changed from UNIT=13 to UNIT=14 to write out 
!                  extcurv.out when doing slab calculations
!                  [unit 13 being used for heating_rate output]
! 20.06.30 (BTD) * correction: added factor 0.921 to calculation of DTAU
! 20.09.10 (BTD) * remove CMB from MMP_RADIATION_FIELD
!                * add CMB to all starlight models
! 20.09.11 (BTD) * replace EXT_F99 with EXT_HD20
! 20.09.20 (BTD) * correct so that U=1 corresponds to mMMP heating rate
! 21.01.07 (BTD) * add option BPASS to use BPASS stellar population models
! 22.12.07 (BTD) * Cosmetic changes
!                * correct value of ICGREMT2 to accomodate DH21c astrodust
!                  dielectric functions (dependent on shape, porosity,
!                  and metallic Fe fraction)
! end history
! --------------------------------------------------------------------

      ICGREMT2=29+540   ! = 569: D16 graphite, MG EMT with matrix E perp c

      PI=4.D0*ATAN(1.D0)
      H=6.62607D-27
      C=2.997925D10
      HC=H*C
      TCMB=2.725

!*** diagnostic
      write(0,*)'vsg_stat_therm_v12 ckpt 0'
!***
      OPEN(UNIT=13,FILE='heating_rate.out')
      WRITE(13,9713)
 9713 FORMAT('  a(cm)   P(erg s-1) comp. b/a Ndot(s-1) <hnu>/eV')

!------------ wavelength grid for illuminating radiation ---------------

      LAMBDAMIN=0.0912*1.D-4
      LAMBDAMAX=8000.*1.D-4

! NISRF         = number of wavelengths used for ISRF
! WL(1-NISRF)   = wavelengths (cm)
! ISRF(1-NISRF) = c*u_lambda (erg cm-3 s-1)

      NISRF=2500
      IF(NISRF.GT.NISRF0)THEN
         WRITE(0,*)'ERROR: NISRF=',NISRF,' > NISRF0=',NISRF0
         STOP
      ENDIF
      DLGLAMBDA=DLOG10(LAMBDAMAX/LAMBDAMIN)/(NISRF-1.)

!-----------------specify whether single zone or multizone ------------------
! BSH: Reducing to just ISRF

      OPEN(UNIT=3,FILE='vsg_zones.par',STATUS='OLD')
      READ(3,*)GEOM

      IF(GEOM.EQ.'isrf')GEOM='ISRF'

      IF(GEOM.EQ.'ISRF')THEN
         WRITE(0,FMT='(A,A)')GEOM,' = GEOM'
         NZONE=1
      ELSE
         WRITE(0,*)'ERROR: unrecognized option GEOM=',GEOM
         STOP
      ENDIF
      CLOSE(3)

      DO I=1,NISRF
         ISRF_WL(I)=LAMBDAMIN*10.**((I-1.)*DLGLAMBDA)
      ENDDO

!--------------------- specify the radiation field ---------------------
! BSH: Removing all but MMP
      
!     radiation='MMPISRF' for MMP ISRF
      
!*** diagnostic
      write(0,fmt='(a)')
     &   'vsg_stat_therm_v12 ckpt 1, about to read vsg_isrf.par'
!***
      OPEN(UNIT=3,FILE='vsg_isrf.par',STATUS='OLD')
      READ(3,*)RADIATION
      IF(RADIATION.EQ.'MMPISRF')THEN

! UMMP = u(starlight)/u(mMMP83 starlight)

         READ(3,*)UMMP
         
! add other radiation field options here:

      ELSE
         WRITE(0,fmt='(a,a,a)')'vsg_stat_therm_v12 ckpt 2 fatal ',
     &      'error: unrecognized radiation type=',RADIATION
         STOP
      ENDIF
      CLOSE(3)

!*** diagnostic
      write(0,fmt='(a)')
     &   'vsg_stat_therm_v12 ckpt 3, completed reading vsg_isrf.par'
!***
! specify grain component:
!	composition='astrosil' (Draine 2003 astrosilicate)
!                   'gra_PAHi' (graphite-PAH+)
!                   'gra_PAHn' (graphite-PAH0)
!                   'gra_PAHs' (graphite-PAH0 and PAH+)
!	            'sil_hd14' (Hensley+Draine 2014 astrosil)
!                   'sil_Fe02' (Hensley+Draine 2014 astrosil + 2% of Fe)
!                   'sil_Fe05' (Hensley+Draine 2014 astrosil + 2% of Fe)
!                   'sil_Fe10' (Hensley+Draine 2014 astrosil + 2% of Fe)
!                   'sil_Fe20' (Hensley+Draine 2014 astrosil + 2% of Fe)
!                   'Fe__dh13' (Draine+Hensley 2013 metallic Fe)
!                   'gPAHi_16' (Draine 2016 graphite-PAH+, at=20A)
!                   'gPAHn_16' (Draine 2016 graphite-PAH0, at=20A)
!                   'gPAHi16b' (Draine 2016 graphite-PAH+, at=50A)
!                   'gPAHi16b' (Draine 2016 graphite-PAH0, at=50A)
!                   'DH20Ad_Px.xx_z.zz_ba' (Draine & Hensley 2020 Astrodust)
!                        porosity Px.xx (x.xx = 0.00, 0.10, 0.20, ..., 0.90)
!                        Fe metal fraction = z.zz (=0.00, 0.10, 0.20, ..., 0.50)
!                        axial ratio ba
!                        (ba = 0.333,0.400,0.500,0.556,0.625,0.714,0.833,1.000,
!                              1.200,1.400,1.600,1.800,2.000,2.500,3.000)
      
! specify number of grain sizes, a_min (cm) , a_max (cm)

!*** diagnostic
      write(0,fmt='(a)')
     &   'vsg_stat_therm_v12 ckpt 4, about to read vsg_size.par'
!***
      OPEN(UNIT=3,FILE='vsg_size.par',STATUS='OLD')
      READ(3,*)COMPOSITION
      READ(3,*)ICASE
      READ(3,*)ISHAP   ! -7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7
      READ(3,*)NSIZE
      IF(NSIZE.GT.NSIZE0)THEN
         WRITE(0,*)'Error: vsg_size.par calls for NSIZE=',NSIZE
         WRITE(0,*)'exceeding maximum NSIZE0=',NSIZE0
         STOP
      ENDIF
      READ(3,*)AMIN
      READ(3,*)AMAX
      CLOSE(3)

!*** diagnostic
      write(0,fmt='(a)')
     &   'vsg_stat_therm_v12 ckpt 5, done reading vsg_size.par'
!***

!*** diagnostic
      write(0,fmt='(a)')
     &   'vsg_stat_therm_v12 ckpt 5.1, about to read vsg_mrf.par'
!***
      OPEN(UNIT=3,FILE='vsg_mrf.par',STATUS='OLD')
      READ(3,*)DELTAMIN
      READ(3,*)DELTAMAX
      READ(3,*)LAMBD_ULAMBD
      CLOSE(3)

!*** diagnostic
      write(0,fmt='(a)')
     &   'vsg_stat_therm_v12 ckpt 5.2, done reading vsg_mrf.par'
!***

!------------ begin sanity check -----------------------------
      PORO=0.
      IF(COMPOSITION.EQ.'astrosil_D2003')THEN
         ICAS=1
      ELSEIF(COMPOSITION.EQ.'gra_PAHi_D2003')THEN
         ICAS=2
      ELSEIF(COMPOSITION.EQ.'gra_PAHn_D2003')THEN
         ICAS=3
      ELSEIF(COMPOSITION.EQ.'sil_HD14_Fe_00')THEN
         ICAS=4
      ELSEIF(COMPOSITION.EQ.'sil_HD14_Fe_02')THEN
         ICAS=5
      ELSEIF(COMPOSITION.EQ.'sil_HD14_Fe_05')THEN
         ICAS=6
      ELSEIF(COMPOSITION.EQ.'sil_HD14_Fe_10')THEN
         ICAS=7
      ELSEIF(COMPOSITION.EQ.'sil_HD14_Fe_20')THEN
         ICAS=8
      ELSEIF(COMPOSITION.EQ.'Fe_metal__DH13')THEN
         ICAS=9
      ELSEIF(COMPOSITION.EQ.'gra_PAHi__2014')THEN
         ICAS=10
      ELSEIF(COMPOSITION.EQ.'gra_PAHn__2014')THEN
         ICAS=11
      ELSEIF(COMPOSITION.EQ.'gra_PAHib_2014')THEN
         ICAS=12
      ELSEIF(COMPOSITION.EQ.'gra_PAHnb_2014')THEN
         ICAS=13
      ELSEIF(COMPOSITION.EQ.'gra_D16rc_PAHi')THEN
         ICAS=14
      ELSEIF(COMPOSITION.EQ.'gra_D16rc_PAHn')THEN
         ICAS=15
      ELSEIF(COMPOSITION.EQ.'gra_D16rcPAHib')THEN
         ICAS=16
      ELSEIF(COMPOSITION.EQ.'gra_D16rcPAHnb')THEN
         ICAS=17
      ELSEIF(COMPOSITION.EQ.'gra_D16ca_PAHi')THEN
         ICAS=18
      ELSEIF(COMPOSITION.EQ.'gra_D16ca_PAHn')THEN
         ICAS=19
      ELSEIF(COMPOSITION.EQ.'gra_D16caPAHib')THEN
         ICAS=20
      ELSEIF(COMPOSITION.EQ.'gra_D16caPAHnb')THEN
         ICAS=21
      ELSEIF(COMPOSITION.EQ.'gra_D16emtPAHi')THEN
         ICAS=22
      ELSEIF(COMPOSITION.EQ.'gra_D16emtPAHn')THEN
         ICAS=23
      ELSEIF(COMPOSITION.EQ.'graD16emtPAHib')THEN
         ICAS=24
      ELSEIF(COMPOSITION.EQ.'graD16emtPAHnb')THEN
         ICAS=25
      ELSEIF(COMPOSITION.EQ.'graphite_D16rc')THEN
         ICAS=26
      ELSEIF(COMPOSITION.EQ.'graphite_D16ca')THEN
         ICAS=27
      ELSEIF(COMPOSITION.EQ.'graphiteD16emt')THEN
         ICAS=28
      ELSEIF(COMPOSITION.EQ.'graphD16emtMG2')THEN
         ICAS=ICGREMT2

! ICASE = 29 + IDX
!              IDX = 6*JPO + JFE = 0 to 90
!                 JPO=0-9 for Poro=0, 0.1, 0.2, ..., 0.9
!                 JFE=0-5 for f_Fe=0, 0.1, 0.2, 0.3, 0.4, 0.5

      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.00_0.00')THEN
         ICAS=29
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.00_0.10')THEN
         ICAS=29+1
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.00_0.20')THEN
         ICAS=29+2
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.00_0.30')THEN
         ICAS=29+3
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.00_0.40')THEN
         ICAS=29+4
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.00_0.50')THEN
         ICAS=29+5

      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.10_0.00')THEN
         ICAS=29+6
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.10_0.10')THEN
         ICAS=29+7
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.10_0.20')THEN
         ICAS=29+8
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.10_0.30')THEN
         ICAS=29+9
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.10_0.40')THEN
         ICAS=29+10
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.10_0.50')THEN
         ICAS=29+11

      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.20_0.00')THEN
         ICAS=29+12
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.20_0.10')THEN
         ICAS=29+13
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.20_0.20')THEN
         ICAS=29+14
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.20_0.30')THEN
         ICAS=29+15
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.20_0.40')THEN
         ICAS=29+16
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.20_0.50')THEN
         ICAS=29+17

      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.30_0.00')THEN
         ICAS=29+18
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.30_0.10')THEN
         ICAS=29+19
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.30_0.20')THEN
         ICAS=29+20
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.30_0.30')THEN
         ICAS=29+21
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.30_0.40')THEN
         ICAS=29+22
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.30_0.50')THEN
         ICAS=29+23

      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.40_0.00')THEN
         ICAS=29+24
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.40_0.10')THEN
         ICAS=29+25
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.40_0.20')THEN
         ICAS=29+26
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.40_0.30')THEN
         ICAS=29+27
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.40_0.40')THEN
         ICAS=29+28
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.40_0.50')THEN
         ICAS=29+29
         
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.50_0.00')THEN
         ICAS=29+30
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.50_0.10')THEN
         ICAS=29+31
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.50_0.20')THEN
         ICAS=29+32
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.50_0.30')THEN
         ICAS=29+33
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.50_0.40')THEN
         ICAS=29+34
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.50_0.50')THEN
         ICAS=29+35
        
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.60_0.00')THEN
         ICAS=29+36
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.60_0.10')THEN
         ICAS=29+37
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.60_0.20')THEN
         ICAS=29+38
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.60_0.30')THEN
         ICAS=29+39
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.60_0.40')THEN
         ICAS=29+40
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.60_0.50')THEN
         ICAS=29+41
         
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.70_0.00')THEN
         ICAS=29+42
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.70_0.10')THEN
         ICAS=29+43
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.70_0.20')THEN
         ICAS=29+44
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.70_0.30')THEN
         ICAS=29+45
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.70_0.40')THEN
         ICAS=29+46
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.70_0.50')THEN
         ICAS=29+47
         
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.80_0.00')THEN
         ICAS=29+48
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.80_0.10')THEN
         ICAS=29+49
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.80_0.20')THEN
         ICAS=29+50
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.80_0.30')THEN
         ICAS=29+51
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.80_0.40')THEN
         ICAS=29+52
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.80_0.50')THEN
         ICAS=29+53

      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.90_0.00')THEN
         ICAS=29+54
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.90_0.10')THEN
         ICAS=29+55
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.90_0.20')THEN
         ICAS=29+56
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.90_0.30')THEN
         ICAS=29+57
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.90_0.40')THEN
         ICAS=29+58
      ELSEIF(COMPOSITION.EQ.'DH21Ad_P0.90_0.50')THEN
         ICAS=29+59
      ELSEIF(COMPOSITION.EQ.'gra_D16emt2_pahib')THEN
         ICAS=90
      ELSEIF(COMPOSITION.EQ.'gra_D16emt2_pahnb')THEN
         ICAS=91

      ELSE
         WRITE(0,fmt='(a,a,a)')'vsg_stat_therm_v12 ckpt 6 fatal error:',
     &      ' unrecognized composition=',COMPOSITION
         STOP
      ENDIF
      IF(ICAS.NE.ICASE)THEN
         WRITE(0,fmt='(a,a,a)')
     &      'vsg_stat_therm_v12 ckpt 7 fatal error: ',
     &      'COMPOSITION=',COMPOSITION
         WRITE(0,*)'not consistent with ICASE=',ICASE
         STOP
      ENDIF

! set porosity

      PORO=0.
      IF(ICAS.GE.29.AND.ICAS.LT.ICGREMT2)PORO=0.10*INT((ICAS-29)/6.)

      IF(ISHAP.EQ.-7)THEN
         BOVERA=0.333
      ELSEIF(ISHAP.EQ.-6)THEN
         BOVERA=0.400
      ELSEIF(ISHAP.EQ.-5)THEN
         BOVERA=0.500
      ELSEIF(ISHAP.EQ.-4)THEN
         BOVERA=0.556
      ELSEIF(ISHAP.EQ.-3)THEN
         BOVERA=0.625
      ELSEIF(ISHAP.EQ.-2)THEN
         BOVERA=0.714
      ELSEIF(ISHAP.EQ.-1)THEN
         BOVERA=0.833
      ELSEIF(ISHAP.EQ.0)THEN
         BOVERA=1.000
      ELSEIF(ISHAP.EQ.1)THEN
         BOVERA=1.200
      ELSEIF(ISHAP.EQ.2)THEN
         BOVERA=1.400
      ELSEIF(ISHAP.EQ.3)THEN
         BOVERA=1.600
      ELSEIF(ISHAP.EQ.4)THEN
         BOVERA=1.800
      ELSEIF(ISHAP.EQ.5)THEN
         BOVERA=2.000
      ELSEIF(ISHAP.EQ.6)THEN
         BOVERA=2.500
      ELSEIF(ISHAP.EQ.7)THEN
         BOVERA=3.000
      ELSE
         WRITE(0,FMT='(A,I5)')
     &      'vsg_stat_therm_v12 ckpt 0: invalid ISHAP=',ISHAP
         STOP
      ENDIF

!------------------------ end sanity check -------------------------

! specify number of energy bins to use in calculation
! specify method = 'stati' for "exact-statistical" method
!                  'dbdis' for "thermal-discrete" method
!                  'dbcon' for "thermal-continous" method

!*** diagnostic
      write(0,fmt='(a)')
     &   'vsg_stat_therm_v12 ckpt 8, about to read vsg_meth.par'
!***
      OPEN(UNIT=3,FILE='vsg_meth.par',STATUS='OLD')
      READ(3,*)METHOD
      READ(3,*)NSTATE
      CLOSE(3)
!*** diagnostic
      write(0,fmt='(a)')
     &   'vsg_stat_therm_v12 ckpt 9, completed reading vsg_meth.par'
      write(0,*)
     &   '   NSTATE=',NSTATE
!***
! Sanity check:

      IF(NSTATE.GT.NSTATEMX)THEN
         WRITE(0,*)'ERROR: called for NSTATE=',NSTATE,
     &             ' > NSTATEMX=',NSTATEMX
         STOP
      ENDIF

! Calculate some ISRF parameters and write them out

      IF(RADIATION.EQ.'MMPISRF')THEN
         WRITE(CDESCR(1),FMT='(A)')' mMMP83 starlight'
         WRITE(CDESCR(2),FMT='()')
      ENDIF

      DO I=1,NISRF

         IF(RADIATION.EQ.'MMPISRF')THEN
            CALL MMP_RADIATION_FIELD(ISRF_WL(I),ISRF(I))
!*** diagnostic
!            if(i.eq.300)write(0,fmt='(a,1pe10.3)')
!     &         'vsg_stat_therm_v12 ckpt 9.5: isrf(300)=',isrf(300)
!***
            ISRF(I)=UMMP*ISRF(I)
            WRITE(RADDESCRIPT,FMT='(A20,F7.3)')
     &         ' MMP83 spec, log(U)=',LOG10(UMMP)
         
         ENDIF

! add CMB:

         CALL BB_RADIATION_FIELD(TCMB,1.D0,ISRF_WL(I),ICMB)
         ISRF(I)=ISRF(I)+ICMB

! if desired, can add other components (e.g., dust emission) here

!*** diagnostic
!         write(0,fmt='(a,1p2e10.3,a)')
!     &      'vsg_stat_therm_v12 ckpt 10:',isrf_wl(i),isrf(i)
!***
      ENDDO   ! enddo i=1,nisrf

      UCMB=7.56577E-15*TCMB**4

!*** diagnostic
         write(0,fmt='(a)')
     &      'vsg_stat_therm_v12 ckpt 10.5:'
!***

! calculate energy density in radiation field

      UISRF=0.
      DO I=1,NISRF
         WGT=1.
         IF(I.EQ.1.OR.I.EQ.NISRF)WGT=0.5
         UISRF=UISRF+WGT*ISRF_WL(I)*ISRF(I)
      ENDDO
      UISRF=UISRF*DLGLAMBDA*LOG(10.)/3.E10
      USTARLIGHT=UISRF-UCMB
      WRITE(0,FMT='(A,A)')RADIATION,' = starlight spectrum'
      WRITE(0,FMT='(1PE10.3,A)')UMMP,' = U (heating rate factor)'
      WRITE(0,FMT='(1PE10.3,A)')USTARLIGHT,
     &                          ' = starlight energy density (erg cm-3)'
      WRITE(0,FMT='(1PE10.3,A)')UCMB,' = CMB energy density (erg cm-3)'
      WRITE(0,FMT='(1PE10.3,A)')UISRF,
     &                          ' = total energy density (erg cm-3)'
      
! write out radiation field to file isrf.out

      OPEN(UNIT=14,FILE='isrf.out')
      WRITE(14,FMT='(A)')RADDESCRIPT
      WRITE(14,FMT='(1PE10.3,A)')UMMP,' = U (heating rate factor)'
      WRITE(14,FMT='(1PE10.3,A)')UISRF,' = u( h nu < 13.6eV) (erg cm-3)'
      WRITE(14,FMT='(1PE10.3,A)')USTARLIGHT,
     &                       ' = u_starlight( h nu < 13.6eV) (erg cm-3)'
      SUM=0.
      SUMD=0.
      DO IWAV=1,NISRF
         SUM=SUM+ISRF_WL(IWAV)*ISRF(IWAV)
         SUMD=SUMD+ISRF_WL(IWAV)**2*ISRF(IWAV)
      ENDDO
      EBAR=1.23984e-4*SUM/SUMD
      WRITE(14,FMT='(1PE10.3,A)')EBAR,' = <h nu> (eV)'
      WRITE(14,FMT='(A,/,A)')
     &   ' lambda lambda*u_lambda',
     &   '   (um)    (erg cm-3)'
      DO IWAV=1,NISRF
         IF(ISRF_WL(IWAV).LE.DELTAMIN.OR.ISRF_WL(IWAV).GE.DELTAMAX)THEN ! put loop in here over lambda and set ISRF to zero
            ISRF(IWAV) = 0.000
            WRITE(14,FMT='(1PE10.3,1PE10.3)')
     &   1.E4*ISRF_WL(IWAV),ISRF_WL(IWAV)*ISRF(IWAV)/2.997925E10
         ELSE
            ISRF(IWAV) = LAMBD_ULAMBD
            WRITE(14,FMT='(1PE10.3,1PE10.3)')
     &   1.E4*ISRF_WL(IWAV),ISRF_WL(IWAV)*ISRF(IWAV)/2.997925E10
         ENDIF
      ENDDO
      CLOSE(14)

!*** end calculation of energy density

      WRITE(CDESCR(6),'(A,A)')METHOD,' = method'
      WRITE(CDESCR(7),'(I5,A)')NSTATE,' = N_bins'
      WRITE(CDESCR(8),'(I5,A)')NISRF,' = N_isrf'
      WRITE(CDESCR(9),'(A)')' '

! ----------------------------------------------------------------

! amin,amax: size in cm;
! NEED TO SPECIFY **nsize** and **nstate**;
! NB: with AMAX=5.01187 um = 10**5.700 A
!          AMIN=3.548   A  = 10**0.550 A
!
!         NSIZE= 84  -> Delta(log10(a))=0.05
!         NSIZE=167                     0.025

      IF(NSIZE.GT.1)THEN
         DELTA_LGA=DLOG10(AMAX/AMIN)/(NSIZE-1.0)
      ELSE
         DELTA_LGA=0.
      ENDIF

      NMAX=NSTATE**2+1

!*** diagnostic
!      write(0,*)'vsg_stat_therm_v12 ckpt 11: NMAX=',NMAX
!      write(0,*)'                           NMAX0=',NMAX0
!      write(0,*)'                 (should have NMAX0.ge.NMAX)'
!*** 

! specify method:
!	method='stati'
!	method='dbdis'
!	method='dbcon'
!	method='dbgdl'

! vsg_stat_therm.pout is output file with energy probability distribution
! vsg_stat_therm.iout is output file with emission spectrum

      OPEN(UNIT=66,FILE='vsg_stat_therm.pout')
      OPEN(UNIT=70,FILE='vsg_stat_therm.dpdtlib.out')
      OPEN(UNIT=77,FILE='vsg_stat_therm.iout')

      IF(GEOM.EQ.'ISRF')THEN
         WRITE(66,FMT='(A,/,A,A,/,A)')
     &      'diffuse starlight with intensity such that',
     &      'aeff=0.1um Astrodust grain (P=0.2, b/a=0.5) absorbs power',
     &      ' h = U*1.829e-12 erg/s',
     &      'P_j = probability of being in energy bin j'
         WRITE(77,FMT='(A,/,A,A,/,A)')
     &      'diffuse starlight with intensity such that',
     &      ' 0.1um Astrodust grain (P=0.2, b/a=0.5) absorbs power',
     &      ' h = U*1.829e-12 erg/s',
     &      'P = time-averaged radiated power/grain'
      ENDIF

! define ISRF0 to save value of unreddened ISRF

      DO IWAV=1,NISRF
         ISRF0(IWAV)=ISRF(IWAV)
      ENDDO

! loop over grain size -----------------------------------------------------

      DO ISIZE=1,NSIZE

         A(ISIZE)=AMIN*10.**((ISIZE-1)*DELTA_LGA)

         write(0,9940)A(ISIZE)
 9940    format(80('='),/,'     a =',1pe10.3,' cm')

         
         write(0,*)'================================================'
         write(0,fmt='(a,a)')'vsg_stat_therm_v12 ckpt 12, ',
     &      ' about to call vsg_td_emission'
         write(0,9731)a(isize),COMPOSITION,UMMP
         write(0,*)'================================================'
         NSTATE1=NSTATE-1

         IF(GEOM.EQ.'ISRF')THEN

!*** diagnostic
            write(0,fmt='(a)')'vsg_stat_therm_v12 ckpt 12.1'
!***
            IWRITE=1
            CALL VSG_TD_EMISSION(IWRITE,METHOD,NCASEMX,ICASE,BOVERA,
     &                           CASEDESCRIPT,RADDESCRIPT,
     &                           DTOH,A(ISIZE),PORO,
     &                           NISRF,ISRF_WL,ISRF,DLGLAMBDA,
     &                           CABS,CSCA,NSTATE,CDESCR,
     &                           U,UA,UB,LNGU,T,DELTA_T,PSTATE,
     &                           RATE_HEATING,COOLING_TIME,RAD_TIME,
     &                           AMATRIX,NSTATE1,PSTATE1,RHS,
     &                           AMATRIX1,NMAX,IJA,SA,
     &                           EMISSION,ABSORPTION,PCUMULT,BGD,
     &                           NDIME,NDIMMODES,
     &                           E,EMODES,LNG,LNGSUM,SCRW)

!*** diagnostic
            write(0,fmt='(a)')'vsg_stat_therm_v12 ckpt 12.2'
!***
         ENDIF

         WRITE(66,6601)(CDESCR(J),J=1,9)
!*** diagnostic
         write(0,fmt='(a)')'vsg_stat_therm_v12 ckpt 20.1'
!***
         WRITE(70,FMT='(I5,A,F9.2,A,A25,A30,A,1PE10.3)')
     &      NSTATE1,' levels:',1.E8*A(ISIZE),' A ',CASEDESCRIPT,
     &      RADDESCRIPT,' Pheat=',RATE_HEATING
!*** diagnostic
         write(0,fmt='(a)')'vsg_stat_therm_v12 ckpt 20.2'
!***
         WRITE(70,FMT='(A,A)')'  T(K) E/hc(cm-1) dE/hc(cm-1) ',
     &                        ' Prob  dE/dT(erg/K)'
!*** diagnostic
         write(0,fmt='(a)')'vsg_stat_therm_v12 ckpt 20.2'
!***
         DO J=1,NSTATE1
            IF(NSTATE1.GT.1)THEN
               IF(J.EQ.1)THEN
                  DEDT=(U(2)-U(1))/(T(2)-T(1))
               ELSEIF(J.EQ.NSTATE1)THEN
                  DEDT=(U(NSTATE1)-U(NSTATE1-1))/
     &                 (T(NSTATE1)-T(NSTATE1-1))
               ELSE
                  DEDT=(U(J+1)-U(J-1))/(T(J+1)-T(J-1))
               ENDIF
            ENDIF
            IF(PSTATE1(J).GE.0.D0)THEN
               IF(T(J).LT.1.E4)THEN
                  WRITE(66,6660)T(J),U(J)/HC,(UB(J)-UA(J))/HC,LNGU(J),
     &                          PSTATE(J),PCUMULT(J),COOLING_TIME(J),
     &                          RAD_TIME(J)
               ELSE
                  WRITE(66,6670)T(J),U(J)/HC,(UB(J)-UA(J))/HC,LNGU(J),
     &                          PSTATE(J),PCUMULT(J),COOLING_TIME(J),
     &                          RAD_TIME(J)
               ENDIF
            ELSE
               IF(T(J).LT.1.E4)THEN
                  WRITE(66,6661)T(J),U(J)/HC,(UB(J)-UA(J))/HC,LNGU(J),
     &                          PSTATE(J),PCUMULT(J),COOLING_TIME(J),
     &                          RAD_TIME(J)
               ELSE
                  WRITE(66,6671)T(J),U(J)/HC,(UB(J)-UA(J))/HC,LNGU(J),
     &                          PSTATE(J),PCUMULT(J),COOLING_TIME(J),
     &                          RAD_TIME(J)
               ENDIF
            ENDIF
            IF(T(J).LT.1.E4)THEN
               WRITE(70,7010)T(J),U(J)/HC,(UB(J)-UA(J))/HC,PSTATE(J),
     &                       DEDT
            ELSE
               WRITE(70,7011)T(J),U(J)/HC,(UB(J)-UA(J))/HC,PSTATE(J),
     &                       DEDT
            ENDIF
         ENDDO   ! enddo j=1,nstate1
      ENDDO   ! enddo isize=1,nsize
      CLOSE(66)
      CLOSE(70)
      CLOSE(77)
      WRITE(0,FMT='(A)')'vsg_stat_therm_v12 ckpt 99: normal termination'
      STOP
 6601 FORMAT(9(A,/),
     &       '     T_j    E_j/hc   deltaE_j   ln(g_j)   P_j',
     &       '     sum(P_j)   tau_cool  tau_rad',/,
     &       '      K      cm-1      cm-1                  ',
     &       '                  sec       sec')
 6660 FORMAT(F9.3,1P7E10.3)
 6661 FORMAT(F9.3,1P3E10.3,1P2E10.2,1P2E10.3)
 6670 FORMAT(F9.2,1P7E10.3)
 6671 FORMAT(F9.2,1P3E10.3,1P2E10.2,1P2E10.3)
 7010 FORMAT(F8.3,1P4E10.3)
 7011 FORMAT(F8.2,1P4E10.3)
 9731 format(50('='),/,5x,
     &       'find dp/dT for a=',1pe10.3,' cm ',A,/,
     &       5X,' (htg rate)/(MMP83 htg rate)=',
     &       1PE10.3,' [for aeff=0.1um DH21Ad_P0.20_0.00_0.500]'/,
     &       50('='))
      END
