      SUBROUTINE SIL_GRA_PAH_TEMP(ICASE,BOVERA,A,NC,NH,ND,
     &                            NAT,DFRAC,NISRF,ISRF_WL,
     &                            ISRF,DLGLAMBDA,TD,
     &                            ED,RATE_HEATING,CABS,CSCA,SCRW,DEDT)
      IMPLICIT NONE
!----------------------- sil_gra_pah_temp_v10 -------------------------
! arguments

      INTEGER ICASE,NISRF
      DOUBLE PRECISION A,BOVERA,DLGLAMBDA,DEDT,DFRAC,ED,
     &   NAT,NC,ND,NH,RATE_HEATING,TD
      DOUBLE PRECISION 
     &   CABS(NISRF),
     &   CSCA(NISRF),
     &   ISRF_WL(NISRF),
     &   ISRF(NISRF),
     &   SCRW(NISRF)

! local variables

      INTEGER I,ICGREMT2,NSTEP
      DOUBLE PRECISION 
     &   CK,EKT,FHI,FLO,RATE_COOLING,RATE_PHOTONS,THI,TLO,TOL

! external functions

      DOUBLE PRECISION AGLI_FZERO,AGLI_ZBRENT_D,PLANCK
      EXTERNAL AGLI_FZERO,AGLI_ZBRENT_D,PLANCK

! ---------------------------------------------------------------------
! subroutine SIL_GRA_PAH_TEMP v10
! given:
!     ICASE = 1 : Draine 2003 astrosilicate
!             2 : graphite-PAH+
!             3 : graphite-PAH0
!             4 : HD14 astrosilicate
!             5 : HD14 astrosilicate with Fe inclusions, f_fill=0.02 
!             6 : HD14 astrosilicate with Fe inclusions, f_fill=0.05 
!             7 : HD14 astrosilicate with Fe inclusions, f_fill=0.10 
!             8 : HD14 astrosilicate with Fe inclusions, f_fill=0.20
!             9 : DH13 metallic Fe
!            10 : graphite-PAH+ as per Hensley & Draine 2014
!                 (reduced sigma_dc(E||c) and reduced a_t=0.0020um)
!            11 : graphite-PAHn as per Hensley & Draine 2014
!                 (reduced sigma_dc(E||c) and reduced a_t=0.0020um)
!            12 : graphite-PAH+ as per Hensley & Draine 2014
!                 (reduced sigma_dc(E||c) and a_t=0.0050um)
!            13 : graphite-PAHn as per Hensley & Draine 2014
!                 (reduced sigma_dc(E||c) and a_t=0.0050um)
!            14 : Draine 2016 graphite-PAH+, random c, a_t=0.0020um 
!            15 : Draine 2016 graphite-PAH0, random c, a_t=0.0020um 
!            16 : Draine 2016 graphite-PAH+, random c, a_t=0.0050um 
!            17 : Draine 2016 graphite-PAH0, random c, a_t=0.0050um 
!            18 : Draine 2016 graphite-PAH+, c || a, a_t=0.0020um 
!            19 : Draine 2016 graphite-PAH0, c || a, a_t=0.0020um 
!            20 : Draine 2016 graphite-PAH+, c || a, a_t=0.0050um 
!            21 : Draine 2016 graphite-PAH0, c || a, a_t=0.0050um 
!            22 : Draine 2016 graphite-PAH+, EMT, a_t=0.0020um 
!            23 : Draine 2016 graphite-PAH0, EMT, a_t=0.0020um 
!            24 : Draine 2016 graphite-PAH+, EMT, a_t=0.0050um 
!            25 : Draine 2016 graphite-PAH0, EMT, a_t=0.0050um 
!            26 : Draine 2016 graphite, random c
!            27 : Draine 2016 graphite, c || a
!            28 : Draine 2016 graphite, EMT (MG1: MG with E||c as matrix)
!       29 -568 : DH21Ad Astrodust (Draine & Hensley 2021a, ApJ 909:94)
!                 ICASE=29+IDX   (IDX from 0 to 59)
!                 IDX = 6*JPO + JFE 
!                     JPO=0-9 for Poro = 0.0, 0.1, ..., 0.9
!                     JFE=0-5 for f_Fe = 0.0, 0.1, ..., 0.5
!           569 : Draine 2016 graphite, EMT (MG2: MG with Eperpc as matrix)
!
!     BOVERA= b/a for spheroid 
!             (= 0.333, 0.400, 0.500, 0.556, 0.625, 0.714, 0.833, 1.000
!                1.200, 1.400, 1.600, 1.800, 2.000, 2.500, 3.000)
!     A     = radius (cm)
!     NC    = number of C atoms [only used if ICASE=2 or 3]
!     ND    =           D        "
!     NH    =           H        "  (does not include D)
!     NAT   = number of atoms [only used if composition = 'sil' or 'Fem'] ??
!     DFRAC = D/H ratio = ^2H/(^1H+^2H)
!             [only used if composition = 'pah']
!     NISRF = number of wavelengths in tabulated ISRF
!     ISRF_WL[] = wavelengths (cm) for tabulated ISRF
!                 (assumed to be uniformly spaced in lg(lambda)
!     ISRF[]    = c*u_lambda (erg cm-3 s-1) for ISRF
!     DLGLAMBDA = increment in log_10(lambda) characterizing ISRF_WL
!     SCRW[]    = scratch space

! returns:

!     CABS[] = C_abs (cm2) at wavelengths ISRF_WL
!     CSCA[] = C_sca (cm2) at wavelengths ISRF_WL
!     TD     = grain temperature (K) for which heating = <cooling>
!     ED     = grain enthalpy (erg) for which heating = <cooling>
!     RATE_HEATING = grain heating rate (erg/s)
!     DEDT   = dED/dT (erg/K)

! Originally written by Aigen Li, Princeton University
! History
! 99.10.11 (AL)  first written
! 00.07.19 (AL)  modified
! 00.10.31 (AL)  modified to calculate only Teq and Eeq
! 00.11.12 (BTD) cosmetic changes, comments added
! 00.11.28 (BTD) added CABS,CSCA,SCRW to argument list
!                add SCRW to argument list for HEATING_RATE and
!                COOLING_RATE
! 06.06.13 (BTD) changed DTOH to DFRAC to avoid confusion
! 07.12.10 (BTD) changed ZG type to INTEGER
! 07.12.11 (BTD) eliminated unneeded arguments A,COMPOSITION
!                in call to HEATING_RATE
!                added support for options 'gra' and 'amc'
!                which assume PAH-like heat capacity
!                eliminated unneeded arguments A,COMPOSITION
!                in call to COOLING_RATE
! 09.08.22 (BTD) modified to allow arbitrarily high values of TEQ
! 11.02.19 (BTD) v3
!                * modified to support composition 'Fem'
! 14.02.19 (BTD) v5
!                * modified to use ICASE instead of COMPOSITION
!                  to specify grain composition
! 14.07.23 (BTD) v6
!                * modified to support 
!                  ICASE=10 : graphite-PAH+ as per HD14 a_t=0.0020um
!                  ICASE=11 : graphite-PAH0 as per HD14 a_t=0.0020um
! 14.07.24 (BTD) * add support for
!                  ICASE=12 : graphite-PAH+ as per HD14 a_t=0.0050um
!                  ICASE=13 : graphite-PAH0 as per HD14 a_t=0.0050um
! 14.09.03 (BTD) v7
!                * add support for
!                  ICASE=14
!                  ICASE=15
!                  ICASE=16
!                  ICASE=17
!                  ICASE=18
! 16.02.08 (BTD) * add support for
!                  ICASE=19, 20, 21, 22
! 16.02.13 (BTD) * add support for
!                  ICASE=23, 24, 25, 26, 27, 28, 29, 30
! 16.02.19 (BTD) * add support for
!                  ICASE=31, 32, 33
! 16.04.26 (BTD) * add argument RATE_PHOTONS to call to subroutine
!                  HEATING_RATE
!                  [now need to use heating_rate_v2.f]
!                * add RATE_PHOTONS to output written to unit 13
! 17.05.03 (BTD) * add support for ICASE=34,35,36,37,38
! 17.10.10 (BTD) * eliminate ISHAP
!                * introduce BOVERA to argument list
!                * add support for ICASE=39,40
! 17.12.25 (BTD) * add support for ICASE=41,42
! 18.03.01 (BTD) v8 modify for new ICASE scheme,
!                   with DH18 silicates from ICASE=29-308
! 18.08.13 (BTD) * modify for new ICASE scheme for
!                  DH18 silicates from ICASE=29-568
! 19.01.21 (BTD) v9 add support for ICASE=569
! 19.08.01 (BTD) v10 modify to support DH19 Astrodust
! 22.12.07 (BTD) * cosmetic changes
!                * correct ICGREMT2 to accomodate DH21a astrodust cases
! 22.12.09 (BTD) v11 modify
!                * add DEDT to argument list
!     end history
!-----------------------------------------------------------------------
      ICGREMT2=29+540  ! 569: D16 graphite MG EMT with matrix = E perp c
      
!*** diagnostic
!      write(0,fmt='(a,1pe10.3,a,i3,a,0pf6.3)')
!     &   'sil_gra_pah_temp_v10 ckpt 0: a=',a,' icase=',icase,
!     &   ' bovera=',bovera
!***

! calculate Cabs (cm^2) [absorption cross section] for lambda[];
! a: cm; cross_section: cm^2;

!*** diagnostic
!      write(0,*)'sil_gra_pah_temp_v10 ckpt 1: about to begin calls to',
!     &          ' SIL_GRA_PAH_CRS_SCT for icase,ishap=',icase,ishap
!      write(0,*)'isrf_wl from',isrf_wl(1),' to',isrf_wl(nisrf)
!***
      DO I=1,NISRF

! irem/src/sil_gra_pah_cs_v5.f uses opt/src/readqlib_v16.f 
! to read precalculated cross sections

         CALL SIL_GRA_PAH_CRS_SCT(ICASE,BOVERA,A,DFRAC,
     &                            ISRF_WL(I),
     &                            CABS(I),CSCA(I))
      ENDDO

!*** diagnostic
      write(0,*)'sil_gra_pah_temp_v10 ckpt 2'
      write(0,*)'completed calls to SIL_GRA_PAH_CRS_SCT'
      write(0,*)'a=',a
      write(0,*)'about to call HEATING_RATE with NISRF=',NISRF
!***

! calculate heating rate (erg s-1) and photoabsorption rate (s-1)

!*** diagnostic
!      write(0,*)'sil_gra_pah_temp_v10 ckpt 2.5, nisrf=',nisrf
!      write(0,fmt='(a,1pe10.3,a,e10.3)')'     isrf_wl(300)=',
!     &   isrf_wl(300),' isrf(300)=',isrf(300)
!***
      CALL HEATING_RATE(NISRF,ISRF_WL,ISRF,DLGLAMBDA,
     &                  CABS,RATE_HEATING,RATE_PHOTONS,SCRW)

!*** diagnostic
!      WRITE(0,*)'sil_gra_pah_temp_v10 ckpt 3:'
!      WRITE(0,FMT='(A,1PE10.3)')
!     &   'returned from HEATING RATE: RATE_HEATING=',rate_heating
!      write(0,*)'a=',a
!***
      WRITE(13,7130)A,RATE_HEATING,ICASE,BOVERA,RATE_PHOTONS,
     &   (RATE_HEATING/(1.60218E-12*RATE_PHOTONS))
 7130 format(1P2E10.3,1X,I3,0PF6.3,1PE10.3,0PF7.3)
      THI=2000.

! proceed to bracket equilibrium temperature by binary chop

      NSTEP=1
 10   CONTINUE

!*** diagnostic
!      write(0,*)'sil_gra_pah_temp_v10 ckpt 4:'
!      write(0,*)'about to call cooling rate with THI=',THI
!***

      CALL COOLING_RATE(NISRF,ISRF_WL,
     &                  DLGLAMBDA,CABS,THI,RATE_COOLING,SCRW)

!*** diagnostic
!      write(0,*)'sil_gra_pah_temp_v10 ckpt 5'
!      write(0,*)'   returned from cooling rate with rate_cooling=',
!     &          rate_cooling
!***

      FHI=(RATE_COOLING-RATE_HEATING)/(RATE_COOLING+RATE_HEATING)
!*** diagnostic
!      write(0,fmt='(a,1pe10.3)')'sil_gra_pah_temp_v10 ckpt 5.1, fhi=',fhi
!***

      IF(FHI.LT.0.)THEN
         THI=1.2*THI
         GOTO 10
      ENDIF
      TLO=THI/2.

!*** diagnostic
!      write(0,*)'sil_gra_pah_temp_v10 ckpt 6:'
!      write(0,*)'call cooling rate again with TLO=',TLO
!***

      CALL COOLING_RATE(NISRF,ISRF_WL,
     &                  DLGLAMBDA,CABS,TLO,RATE_COOLING,SCRW)

!*** diagnostic
!      write(0,*)'sil_gra_pah_temp_v10 ckpt 7:'
!      write(0,fmt='(a,1pe10.3)')
!     &   'returned from cooling rate with rate_cooling=',rate_cooling
!***

      FLO=(RATE_COOLING-RATE_HEATING)/(RATE_COOLING+RATE_HEATING)

!*** diagnostic
!      write(0,fmt='(a,1pe10.3)')'sil_gra_pah_temp_v10 ckpt 7.1, flo=',flo
!***
      IF((FLO*FHI).GT.0.)THEN
         NSTEP=NSTEP+1
         THI=THI/2.
         GOTO 10
      ELSEIF((FLO*FHI).LE.0.)THEN
         GOTO 20
      ENDIF

 20   TOL=1.D-14

! equilibrium temperature must be in [TLO,THI].
! call zbrent to obtain refined estimate

!*** diagnostic
!      write(0,*)'sil_gra_pah_temp_v10 ckpt 8'
!      write(0,*)'about to call agli_zbrent_d for agli_fzero'
!      write(0,fmt='(a,1p2e10.3)')'with tlo,thi=',tlo,thi
!***
! must now use agli_zbrent_d_v2

      TD=AGLI_ZBRENT_D(A,NISRF,ISRF_WL,
     &                 DLGLAMBDA,CABS,RATE_HEATING,AGLI_FZERO,
     &                 TLO,THI,TOL,SCRW)

!*** diagnostic
!      write(0,*)'sil_gra_pah_temp_v10 ckpt 9'
!      write(0,*)'returned from agli_zbrent_d with TD=',TD
!***
      IF(ICASE.EQ.2.OR.
     &   ICASE.EQ.3.OR.
     &   (ICASE.GE.10.AND.ICASE.LE.28).OR.
     &   ICASE.EQ.ICGREMT2)THEN

!*** diagnostic
!         write(0,*)'sil_gra_pah_temp_v10 ckpt 10:',
!     &   ' about to call PAH_SPEC_HEAT...'
!***
         CALL PAH_SPEC_HEAT(1,NC,NH,ND,TD,EKT,CK)
!*** diagnostic
!         write(0,*)'returned from PAH_SPEC_HEAT'
!***
      ELSEIF(ICASE.EQ.1.OR.
     &       (ICASE.GE.4.AND.ICASE.LE.8).OR.
     &       (ICASE.GE.29.AND.ICASE.LT.ICGREMT2))THEN
         CALL SIL_SPEC_HEAT(1,NAT,TD,EKT,CK)  ! MODE=1: use Debye-like model
      ELSEIF(ICASE.EQ.9)THEN
         CALL FE_SPEC_HEAT(1,NAT,TD,EKT,CK)   ! MODE=1: use Debye-like model
      ELSE
         WRITE(0,*)'sil_gra_pah_temp_v10: Fatal error'
         WRITE(0,*)' unrecognized icase=',icase
         STOP
      ENDIF

      ED=EKT*1.38065D-16*TD
      DEDT=CK*1.38065D-16

      WRITE(0,7700)a*1.D4,Td,(ED/1.60218D-12)

!*** diagnostic
!      write(0,*)'returning from sil_gra_pah_temp_v10'
!***
      RETURN
 7700 FORMAT('size=',F6.4,' micron: Teq=',F7.2,' E=',1PE10.3,' eV')
      END
