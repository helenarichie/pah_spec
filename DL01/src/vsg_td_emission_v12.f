      SUBROUTINE VSG_TD_EMISSION(IWRITE,METHOD,NCASEMX,ICASE,BOVERA,
     &                           CASEDESCRIPT,RADDESCRIPT,
     &                           DFRAC,AGR,PORO,
     &                           NISRF,ISRF_WL,ISRF,DLGLAMBDA,
     &                           CABS,CSCA,NSTATE,CDESCR,
     &                           U,UA,UB,LNGU,T,DELTA_T,PSTATE,
     &                           RATE_HEATING,COOLING_TIME,RAD_TIME,
     &                           AMATRIX,NSTATE1,PSTATE1,RHS,
     &                           AMATRIX1,NMAX,IJA,SA,
     &                           EMISSION,ABSORPTION,PCUMULT,BGD,
     &                           NDIME,NDIMMODES,
     &                           E,EMODES,LNG,LNGSUM,SCRW)

!------------------------- vsg_td_emission_v11 -----------------------

      IMPLICIT NONE

! Arguments:

      INTEGER NCASEMX
      CHARACTER 
     &   CASEDESCRIPT*25,
     &   CDESCR(NCASEMX)*79,
     &   METHOD*5,
     &   RADDESCRIPT*30
      INTEGER ICASE,IWRITE,NDIME,NDIMMODES,NISRF,NMAX,NSTATE,NSTATE1
      INTEGER IJA(NMAX)
      DOUBLE PRECISION AGR,BOVERA,DEDT,DLGLAMBDA,DFRAC,PORO,RATE_HEATING
      DOUBLE PRECISION 
     &   ABSORPTION(NISRF),
     &   AMATRIX(NSTATE,NSTATE),
     &   AMATRIX1(NSTATE1,NSTATE1),
     &   BGD(NSTATE,NSTATE),
     &   CABS(NISRF),
     &   CSCA(NISRF),
     &   COOLING_TIME(NSTATE),
     &   DELTA_T(NSTATE),
     &   E(0:NDIME),
     &   EMISSION(NISRF),
     &   EMODES(NDIMMODES),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF),
     &   LNG(0:NDIME),
     &   LNGSUM(0:NDIME),
     &   LNGU(NSTATE),
     &   PCUMULT(NSTATE),
     &   PSTATE(NSTATE),
     &   PSTATE1(NSTATE1),
     &   RAD_TIME(NSTATE),
     &   RHS(NSTATE1),
     &   SA(NMAX),
     &   SCRW(NISRF),
     &   T(NSTATE),
     &   U(NSTATE),
     &   UA(NSTATE),
     &   UB(NSTATE)

! Local variables:

      LOGICAL
     &   REFINE
      INTEGER 
     &   F,I,ICGREMT2,IDX,IT_END,ITER,ITMAX,ITOL,
     &   J,JCUT,JPO,JEMAX0,JFE,JPMAX,K,
     &   NATOM_INT,NC_INT,ND_INT,NH_INT,NMODES,NSET

      DOUBLE PRECISION
     &   A1,A2,AESVS,APAH, 
     &   C,CV,
     &   DE,DELTA_E,DJ,
     &   E_ABS,E_EM,EBS,EEQ,EEQSS,EMSN_INTRABIN,EMSN_INTERBIN,ERR,
     &   FAC,FFE,FVSIL,FVFE,GKJ,
     &   H,HC,HCC,HTOC,
     &   KB,NATOM,NC,ND,NH,NHTOT,
     &   PI,PMAX,PMIN_LO,PMIN_UP,PRAD,
     &   SUM_PSTATE,
     &   TEQ,TERM,THRESH,TOL,
     &   UMAX,UMAXHI,UMAXLO,UMAXMIN,UMIN,UMINHI,UMINLO

! --------------------------------------------------------------------
! Subroutine VSG_TD_EMISSION

! NB: in recent years, this routine is only used with METHOD='dbdis'
!     and may no longer be correct for other cases

! input: 
!       METHOD = 'stati'   (for "exact-statistical" treatment)
!                'dbdis'   (for "Debye-thermal-discrete" treatment)
!                'stdst'   (for "steady-state" single-T treatment)
!                           [should also set NSTATE =1 and
!                                            NSTATE1=0]
!       NCASEMX = dimensioning parameter for max number of cases
!       ICASE = 1 for D03 astrosilicate
!               2     graphite-PAH+
!               3     graphite-PAH0
!               4     HD2014 silicate, no Fe inclusion
!               5     HD2014 silicate, f_vol= 2% Fe inclusions
!               6     HD2014 silicate, f_vol= 5% Fe inclusions
!               7     HD2014 silicate, f_vol=10% Fe inclusions
!               8     HD2014 silicate, f_vol=20% Fe inclusions
!               9     DH2013 metallic Fe
!              10     graphite-PAH+ as per Hensley & Draine 2014
!                     reduced sigma_dc(E||c) and reduced a_t=0.0020um
!              11     graphite-PAH0 as per Hensley & Draine 2014
!                     reduced sigma_dc(E||c) and reduced a_t=0.0020um
!              12     graphite-PAH+ as per Hensley & Draine 2014
!                     reduced sigma_dc(E||c) and original a_t=0.0050um
!              13     graphite-PAH0 as per Hensley & Draine 2014
!                     reduced sigma_dc(E||c) and original a_t=0.0050um
!              14     D16 graphite-PAH+, random c axis, a_t=20A
!              15     D16 graphite-PAH0, random c axis, a_t=20A
!              16     D16 graphite-PAH+, random c axis, a_t=50A
!              17     D16 graphite-PAH0, random c axis, a_t=50A
!              18     D16 graphite-PAH+, c || a, a_t=20A
!              19     D16 graphite-PAH0, c || a, a_t=20A
!              20     D16 graphite-PAH+, c || a, a_t=50A
!              21     D16 graphite-PAH0, c || a, a_t=50A
!              22     D16 graphite-PAH+, EMT, a_t=20A
!              23     D16 graphite-PAH0, EMT, a_t=20A
!              24     D16 graphite-PAH+, EMT, a_t=50A
!              25     D16 graphite-PAH0, EMT, a_t=50A
!              26     D16 graphite, random c axis
!              27     D16 graphite, c || a
!              28     D16 graphite, EMT (MG1: MG with Eperpc as inclusion)

!          29-568     DH20 Astrodust
!                         poro from 0 to 0.5:  6 cases, steps of 0.1
!                         f_Fe from 0 to 0.5:  6 cases, steps of 0.1
!                         b/a from 0.33 to 3: 15 cases
!                         total number of cases: 6x6x15=540

!             569     D16 graphite, EMT (MG2: MG with Eperpc as matrix)

!       BOVERA = b/a for spheroid
!                (= 0.333, 0.400, 0.500, 0.556, 0.625, 0.714, 0.833, 1.000,
!                   1.200, 1.400, 1.600, 1.800, 2.000, 2.500, 3.000)
!
!       RADDESCRIPT = A30 string describing assumed radiation field
!       DFRAC   = D/(H+D) ratio (deuterium/total hydrogen ratio) 
!                              (fully deuterated grain has DFRAC=1)
!       AGR    = equal-volume-sphere radius (cm) [volume includes voids]
!       PRORO  = Porosity = fraction of volume contributed by voids
!       NISRF  = number of wavelengths in tabulated ISRF
!       ISRF_WL[1-NISRF] = wavelengths (cm) for tabulated ISRF
!                          (assumed to be uniform in lg(lambda))
!       ISRF[1-NISRF] = c*u_lambda (erg cm-3 s-1) for ISRF
!       NSTATE = number of states (temperature/enthalpy)
!              = number of energy bins
!       DLGLAMBDA = delta(log10(lambda)) for wavelengths ISRF_WL[]

! returns:

!       CDESCR[3,4,5]  = descriptive strings

!       CABS[1-NISRF]  = absorption cross section (cm^2);
!       CSCA[1-NISRF]  = scattering cross section (cm^2) for sil/gra;
!       U[1-NSTATE]    = the mid-bin enthalpy (erg)
!       UA[1-NSTATE]   = lower enthalpy limit for each bin (erg)
!       UB[1-NSTATE]   = upper enthalpy limit for each bin (erg)
!       LNGU[1-NSTATE] = ln(degeneracy) of each bin
!       T[1-NSTATE]    = temperature of each bin (K)
!       DELTA_T[1-NSTATE] = temperature width of each bin (K)
!       PSTATE[1:NSTATE] = probability of finding grain in the 
!                   ith bin, i.e., within [ua, ub];
!                   solution to matrix equation
!                   AMATRIX[f,i]*Pstate[i]=0 for all f
!       COOLING_TIME[1-NSTATE] = U/(dU/dt) (sec) for each bin;
!       RAD_TIME[1-NSTATE] = 1/{sum_over_f Afi (i>f)} (sec);
!       AMATRIX(1:NSTATE,1:NSTATE) = transition matrix (s-1)
!       PCUMULT[1:NSTATE] = cumulative probability distribution;
!       EMISSION[1:NISRF] = F_lambda (erg cm-1 s-1)
!                           where F = power radiated by the grain
!       ABSORPTION[1:NISRF] = C_abs(lambda)*c*u_lambda (erg cm-1 s-1)
!                             at wavelengths ISRF_WL
!       NSTATE1 = NSTATE-1 if returning full distribution PSTATE(1:NSTATE)
!                          after solving for steady-state
!               = 1 if returning delta function distribution
!       PSTATE1[1:NSTATE1]
!       RHS[1:NSTATE1] (right-hand-side);
!       AMATRIX1[1:NSTATE1,1:NSTATE] 
!       where PSTATE1 is the solution to the matrix equation

!       AMATRIX1[1:NSTATE1,1:NSTATE1]*PSTATE1[1:NSTATE1]=RHS[1:NSTATE1]

!       NMAX=(NSTATE1+1)**2+1, 
!       IJA[1:NMAX], SA[1:NMAX]: work arrays used by BiCG method; 
!       ABSORPTION[]: array for calculating absorption rate (erg/s);
!       BGD[1:N,1:N]: working array for the continuous cooling method.
!       SCRW[1:NISRF] = scratch workspace

! PSTATE1(1:NSTATE1) is written (with other data) to unit=66 
!                    (set to vsg_stat_therm.pout by program vsg_stat_therm)
!                    and to unit=70 
!                    (set to vsg-stat_therm.dpdtlib.out by vsg_stat_therm)
! emission spectrum is written to unit=77 
!                    (set to vsg_stat_therm.iout by vsg_stat_therm)
!
!-----------------------------------------------------------------------
! Subroutine VSG_TD_EMISSION calculates the energy (temperature)
! probability distribution function for small grains heated by
! starlight.
! Original code created by Aigen Li, Princeton University
! Subsequent modifications by B.T. Draine, Princeton University
! History:
! 99.09.28 (AL)  First created to compute PAH emission
! 99.10.11 (AL)  Converted to subroutine
! 99.10.18 (AL)  Adopt enthalpy bins
!                include silicates and graphites
! 00.01.27 (AL)  exact-statistical solution option
! 00.03.03 (AL)  add continuous cooling approximation
! 00.03.08 (AL)  allow for D instead of H
! 00.10.02 (AL)  calculate equilibrium temperature (TEQ), 
!                equilibrium enthalpy (EEQ) 
! 00.11.11 (BTD) cosmetic changes; added comments
! 00.11.22 (BTD) rewrote module which iteratively adjusts UMAX
! 00.11.23 (BTD) continue modifying module to adjust UMAX,UMIN
! 00.11.28 (BTD) add E,EMODES,LNG,LNGSUM,SCRW to argument list
!                so that memory allocation can be controlled from
!                main program, and also to reduce redundant
!                allocation of memory
! 05.10.31 (BTD) modified module which iteratively adjusts UMAX
!                so that UMAX always remains at least 13.7 eV
! 05.10.31 (BTD) modified to write output file in format expected
!                by program makespec
! 05.11.01 (BTD) modified handling of UMAX refinement in hopes of
!                increasing speed
! 06.04.15 (BTD) corrected minor error in format (g77 tolerated,
!                ifort objected)
! 06.11.03 (BTD) v2: modified to require UMAX > UMAXMIN
!                with UMAXMIN=13.65 eV
! 06.11.04 (BTD) v2: modified to make condition for terminating
!                iteration on UMAX more stringent
!                modify to print out PSTATE(J=3) in addition to
!                J=1,2,497,498,499
! 07.12.10 (BTD) v2: changed type of ZG to INTEGER
!                modified calculation of number of NC,NH,ND
!                for 'PAH' and 'gra' options
!                added support for 'amc" option
! 07.12.27 (BTD) now stop reducing Umin when we reach 1 cm-1
! 08.02.18 (BTD) modified to use delta-function T distribution
!                when steady-state energy content EEQ > EEQSS
!                with EEQSS/hc currently set to 9e5 cm-1
!                SIL_GRA_PAH_TEMP has been modified to support this
!                by adding RATE_HEATING to the argument list
! 08.03.08 (BTD) changed EEQSS from 9e5cm-1*hc = 111.6 eV
!                              to 100 eV
!                because of failure that occured when trying to solve
!                for dP/dE when E_ss=104.8 eV 
!                (14.13A PAH ion, MMP rad field with U=1e7)
!                added check for possible numerical failure in LINBCG
! 08.03.08 (BTD) changed criterion for setting first try UMIN
!                was: if a < 70A -> UMIN=0
!                     else UMIN=0.2*EEQSS
!                now: if EEQ < 0.1*EEQSS -> UMIN=0
!                     else UMIN=0.2*EEQSS
! 09.08.26 (BTD) added new formatting for output when T(I)>1e4K
! 11.02.23 (BTD) v4: modified to recognize composition = 'fem'
!                    (metallic Fe)
! 11.05.18 (BTD) * add energy conservation check
!                  (halt if fail by more than 3%)
!                * add check on accuracy of solving original equation.
!                  repeat conjugate gradient solver if solution does
!                  not meet accuracy requirement at ckpt 5.9
!                * add easy solution for low-U limit,
!                  with criterion for deciding whether we are in low-U
!                  limit or need to do full solution.
!                * eliminated unused options "dbcon" and "dbgdl'
!                * changed some addresses to make them sequential
! 12.10.05 (BTD) v5
!                * add RADDESCRIPT to argument list
!                * write RADDESCRIPT to each dP/dT result
!                * write heating rate to each dP/dT result
! 14.02.17 (BTD) v6
!                * add ICASE,ISHAP to argument list
!                * add ICASE,ISHAP to argument list of SIL_GRA_PAH_TEMP
!                * add ICASE,ISHAP to argument list of TRANSITION_MATRIX
!                * eliminate variable COMPOSITION
! 14.07.22 (BTD) v6
!                * add support for ICASE=10 and ICASE=11
! 14.07.24 (BTD) * add support for ICASE=12 and ICASE=13
! 14.09.05 (BTD) * add support for ICASE=14, 15, 16, 17, 18
! 14.09.10 (BTD) * correct bug for ICASE=14
! 16.02.08 (BTD) * add support for ICASE=19, 20, 21, 22
! 16.02.15 (BTD) v7
!                * add support for ICASE=23-30
! 16.02.20 (BTD) * add support for ICASE=31-33
! 16.10.03 (BTD) * ensure that NATOM is defined
!                * modify to set EEQSS=120 eV if T_eq < 4K
!                  in order to avoid numerical problems
! 17.05.03 (BTD) * add support for ICASE=34-38
!                * change CASEDESCRIPT*15 -> CASEDESCRIPT*24
! 17.05.15 (BTD) * reduce EEQSS for TEQ < 4.7K to try to avoid
!                  numerical problem that was encountered for
!                  U=1e-3, a=1000A, sil_DH17_9000_0.20_Fe_0.00
! 17.06.04 (BTD) * set
!                     NSTATE1=NSTATE-1 if returning PSTATE(1:NSTATE)
!                     NSTATE1=1        if using delta function PSTATE(1)
!                * cosmetic changes
! 17.12.25 (BTD) * add support for ICASE=41,42
! 18.02.28 (BTD) v8 change to be consistent with new ICASE scheme
!                * adjusted calculation of NATOM for silicates
!                  with Fe inclusions
! 18.08.12 (BTD) v9 add PORO to argument list so as to be able
!                   to calculate heat capacity given radius
!                * change to new scheme for DH18 sil
!                  for ICASE=29 - 568
!                * new variable to support treatment of DH18 sil porosity
!                  AESVS = radius of equal-solid-volume-sphere
! 18.10.27 (BTD) * corrected logical error in handling ICASE<29
! 19.01.21 (BTD) v10 ensure compatibiity with ICASE=569
! 19.08.01 (BTD) v11 new ICASE scheme for DH19 Astrodust
! 20.03.06 (BTD) * update to DH20 Astrodust
! 20.05.08 (BTD) * modify to limit increase in Umax when Umax already large
! 20.05.10 (BTD) * add IWRITE to argument list
!                  use IWRITE to enable/disable writing to unit 77
! 20.06.02 (BTD) * improved treatment for increasing UMAX
! 20.06.16 (BTD) * changed assumed RHO_C from 2.26 gm cm-3 to 2.0 gm cm-3
!                  NC from 472 (a/10A)^3 to 418 (a/10A)^3
! 22.12.07 (BTD) * cosmetic changes
!                * change ICGREMT2 from 119 to 569 to accomodate
!                  the 540 cases of optical properties for DH21a astrodust
!                  (shape, porosity, and metallic Fe fraction)
!                * one more place to change rho_C from 2.26 to 2.0
! 22.12.09 (BTD) v12 modify
!                * add DEDT to argument list
!                * add DEDT to argument list of subroutine SIL_GRA_PAH_TEMP
! end history
! -------------------------------------------------------------------- 
! NMODES 40000 is for a=31 angstrom (nc=9397);
! -------------------------------------------------------------------
! call PAH_MODES.f to calculate the exact vibrational modes
! for the exact-statistical method; for other methods just
! arbitrarily assign 0. to EMODES(1:NDIMMODES). 
! 00.03.26 call dens_states.f to calculate the number of states so as
!          to derive level degeneracies. 
! 00.03.26 set umax (erg) at 4*13.6eV. 

!*** diagnostic
      write(0,'(a)')'vsg_td_emission_v11 ckpt 0:'
!      write(0,'(a,i4)')'   nstate =',nstate
!      write(0,'(a,i4)')'   nstate1=',nstate1
!      write(0,'(a,i10)')'   nmax   =',nmax
!      write(0,'(a,i3,a,i3)')'   icase=',ICASE,' ishap=',ISHAP
!      write(0,'(a,1pe10.3)')'   a=',a
!      write(0,'(a,a)')'   raddescript=',raddescript
!      write(0,*)'   nisrf=',nisrf
!      write(0,fmt='(a,1pe10.3)')'   isrf_wl(300)=',isrf_wl(300)
!      write(0,fmt='(a,1pe10.3)')'      isrf(300)=',isrf(300)
!***
      H=6.62607D-27
      C=2.99793D10
      KB=1.38065D-16
      HC=H*C
      HCC=HC*C
      PI=3.141593D0

      UMAXHI=1.D70*HC
      UMAXLO=0.D0
      UMINHI=1.D70*HC
      UMINLO=0.D0

! compute AESVS = radius of equal-solid-volume-sphere

      AESVS=AGR*(1.-PORO)**(1./3.)

      NSTATE1=NSTATE-1

      WRITE(CDESCR(3),FMT='(1PE10.3,A)')AGR,' = a(cm)'
      
      ICGREMT2=29+540 ! 569=case number for D16 graphite, MG EMT (matrix:Eperpc)
      
      IF(ICASE.EQ.1)CASEDESCRIPT='astrosil_D2003'  ! astrosil
      IF(ICASE.EQ.2)CASEDESCRIPT='gra_PAHi_D2003'  ! gra_PAHi
      IF(ICASE.EQ.3)CASEDESCRIPT='gra_PAHn_D2003'  ! gra_PAHn
      IF(ICASE.EQ.4)CASEDESCRIPT='sil_HD14_Fe_00'  ! sil_hd14
      IF(ICASE.EQ.5)CASEDESCRIPT='sil_HD14_Fe_02'  ! sil_Fe02
      IF(ICASE.EQ.6)CASEDESCRIPT='sil_HD14_Fe_05'  ! sil_Fe05 
      IF(ICASE.EQ.7)CASEDESCRIPT='sil_HD14_Fe_10'  ! sil_Fe10
      IF(ICASE.EQ.8)CASEDESCRIPT='sil_HD14_Fe_20'  ! sil_Fe20
      IF(ICASE.EQ.9)CASEDESCRIPT='Fe_metal__DH13'  ! Fe__dh13
      IF(ICASE.EQ.10)CASEDESCRIPT='gra_PAHi__2014' ! gPAHi_14
      IF(ICASE.EQ.11)CASEDESCRIPT='gra_PAHn__2014' ! gPAHn_14
      IF(ICASE.EQ.12)CASEDESCRIPT='gra_PAHib_2014' ! gPAHi14b
      IF(ICASE.EQ.13)CASEDESCRIPT='gra_PAHnb_2014' ! gPAHn14b
      IF(ICASE.EQ.14)CASEDESCRIPT='gra_D16rc_PAHi'
      IF(ICASE.EQ.15)CASEDESCRIPT='gra_D16rc_PAHn'
      IF(ICASE.EQ.16)CASEDESCRIPT='gra_D16rcPAHib'
      IF(ICASE.EQ.17)CASEDESCRIPT='gra_D16rcPAHnb'
      IF(ICASE.EQ.18)CASEDESCRIPT='gra_D16ca_PAHi'
      IF(ICASE.EQ.19)CASEDESCRIPT='gra_D16ca_PAHn'
      IF(ICASE.EQ.20)CASEDESCRIPT='gra_D16caPAHib'
      IF(ICASE.EQ.21)CASEDESCRIPT='gra_D16caPAHnb'
      IF(ICASE.EQ.22)CASEDESCRIPT='gra_D16emtPAHi'
      IF(ICASE.EQ.23)CASEDESCRIPT='gra_D16emtPAHn'
      IF(ICASE.EQ.24)CASEDESCRIPT='graD16emtPAHib'
      IF(ICASE.EQ.25)CASEDESCRIPT='graD16emtPAHnb'
      IF(ICASE.EQ.26)CASEDESCRIPT='graphite_D16rc'
      IF(ICASE.EQ.27)CASEDESCRIPT='graphite_D16ca'
      IF(ICASE.EQ.28)CASEDESCRIPT='graphiteD16emt'
      IF(ICASE.GE.29.AND.ICASE.LT.ICGREMT2)THEN
         IDX=ICASE-29     ! IDX=0,...,539
         JPO=MOD(IDX,6)   ! JPO=0,1,2,3,4,5
         JFE=IDX-6*JPO
         PORO=0.1*JPO
         FFE=0.1*JFE      ! fraction of Fe in metallic form
         WRITE(CASEDESCRIPT,FMT='(A,F4.2,A,F4.2)')
     &      'DH21Ad_P',PORO,'_',FFE
      ENDIF
      IF(ICASE.EQ.ICGREMT2)CASEDESCRIPT='graphD16MG2emt'
      IF(ICASE.EQ.90)CASEDESCRIPT='graD16emt_pahi'
      IF(ICASE.EQ.91)CASEDESCRIPT='graD16emt_pahn'
      IF(ICASE.LT.0.OR.ICASE.GT.ICGREMT2)THEN
         WRITE(0,FMT='(A,I6)')
     &      'vsg_td_emission_v11 fatal error: unknown ICASE=',ICASE
         STOP
      ENDIF

      IF(ICASE.EQ.1.OR.
     &   (ICASE.GE.4.AND.ICASE.LE.8))THEN

! silicate, density 3.5 g cm-3, nominal composition FeMgSiO4
!                               7 atoms, 56+24+28+64=172 amu

         FFE=0.
         IF(ICASE.EQ.5)FFE=0.02
         IF(ICASE.EQ.6.OR.ICASE.EQ.15)FFE=0.05
         IF(ICASE.EQ.7.OR.ICASE.EQ.16)FFE=0.10
         IF(ICASE.EQ.17)FFE=0.15
         IF(ICASE.EQ.8.OR.ICASE.EQ.18)FFE=0.20

         NATOM=(4.*PI/3.)*AESVS**3*(
     &         (1.-FFE)*3.5*7./(172.*1.66D-24)+
     &         FFE*7.87/(55.9*1.66D-24))

      ELSEIF(ICASE.EQ.2.OR.
     &       ICASE.EQ.3.OR.
     &       (ICASE.GE.10.AND.ICASE.LE.13).OR.
     &       ICASE.EQ.90.OR.ICASE.EQ.91)THEN

!*** diagnostic
!         write(0,*)'vsg_td_emission_v11 ckpt 1, a=',a
!*** 

! old code:
! we assume graphite density 2.26 g cm-3
!           carbon density 2.26/(12*1.67e-24)=1.13e23cm-3
!         NC=472.E21*AESVS**3

! 20.06.16 (BTD) change in assumed C density for graphite:
!                now assume 2.0 g cm-3

         NC=418.E21*AESVS**3

         ND=0.
         NH=0.
         IF(NC.LT.1.D4)THEN
            WRITE(CDESCR(4),FMT='(I5,A)')NINT(NC),
     &                                   ' = NC for graphite'
         ELSE
            WRITE(CDESCR(4),FMT='(1PE10.3,A)')NC,
     &         ' = NC for graphite'
         ENDIF

! we assume PAH-like up to radius A1
! graphite-like beyond radius A1
! graphite-like region is assumed to have H:C=0.25

         IF(ICASE.EQ.2.OR.
     &      ICASE.EQ.3.OR.
     &      ICASE.EQ.12.OR.
     &      ICASE.EQ.13.OR.
     &      ICASE.EQ.90.OR.
     &      ICASE.EQ.91)A1=50.D-8
         IF(ICASE.EQ.10.OR.
     &      ICASE.EQ.11)A1=20.D-8
         IF(AGR.LE.A1)THEN
            APAH=AGR
!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 2'
!***
            CALL PAH_SIZE_ATOMS(APAH,NC,NHTOT,HTOC)
!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 3'
!***
            ND=DFRAC*NHTOT
            NH=NHTOT-ND
            IF(NC.LT.1.D4)THEN
               WRITE(CDESCR(4),FMT='(3I5,A)')
     &            NINT(NC),NINT(NH),NINT(ND),' = NC,NH,ND for PAH'
            ELSE
               WRITE(CDESCR(4),FMT='(1P3E10.3,A)')
     &            NC,NH,ND,' = NC,NH,ND for PAH'
            ENDIF
         ELSE
            APAH=A1
!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 4'
!***
            CALL PAH_SIZE_ATOMS(APAH,NC,NHTOT,HTOC)
!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 5'
!***
! we now assume graphite density 2.00 g cm-3
! n_C = 2.00/(12.*1.67e-24) = 9.98e22 cm-3
! N_C = 4.18e23 (a/cm)^3

            ND=DFRAC*NHTOT
            NH=NHTOT-ND
            NC=NC+418.D21*(AESVS**3-A1**3)
            NH=NC+0.25*418.D21*(AESVS**3-A1**3)
         ENDIF
         NATOM=NC+NH+ND

      ELSEIF(ICASE.EQ.9)THEN

! metallic Fe
! density 7.87 g cm-3
! mw = 55.845
! number density 7.87/(55.845*1.6605e-24) = 8.487e22 cm-3
! n^{-1/3} = 2.276 Angstrom
! NATOM = 0.355*(A/1.d-8)^3 = 355*(A/1.d-7)^3

         NATOM=(4.*PI/3.)*AESVS**3*7.87/(55.845*1.6605D-24)

      ELSEIF((ICASE.GE.14.AND.ICASE.LE.28).OR.
     &       ICASE.EQ.569)THEN

! We assume graphite density 2.00 g cm-3
! Up to radius a1 grain is PAH-like

         IF(ICASE.EQ.14.OR.
     &      ICASE.EQ.15.OR.
     &      ICASE.EQ.18.OR.
     &      ICASE.EQ.19.OR.
     &      ICASE.EQ.22.OR.
     &      ICASE.EQ.23)A1=20.D-8
         IF(ICASE.EQ.16.OR.
     &      ICASE.EQ.17.OR.
     &      ICASE.EQ.20.OR.
     &      ICASE.EQ.21.OR.
     &      ICASE.EQ.24.OR.
     &      ICASE.EQ.25)A1=50.D-8
         IF((ICASE.GE.26.AND.ICASE.LE.28).OR.ICASE.EQ.ICGREMT2)A1=0.

         IF(AESVS.LE.A1)THEN
            APAH=AESVS
            CALLPAH_SIZE_ATOMS(APAH,NC,NHTOT,HTOC)
            ND=DFRAC*NHTOT
            NH=NHTOT-ND
            IF(NC.LT.1.D4)THEN
               WRITE(CDESCR(4),FMT='(3I5,A)')
     &            NINT(NC),NINT(NH),NINT(ND),' = NC,NH,ND for PAH'
            ELSE
               WRITE(CDESCR(4),FMT='(1P3E10.3,A)')
     &            NC,NH,ND,' = NC,NH,ND for PAH'
            ENDIF
         ELSE
            APAH=A1
            ND=0.
            NH=0.
            NC=418.D21*AESVS**3
            IF(APAH.GT.0.)THEN
!*** diagnostic
!               write(0,*)'vsg_td_emission_v11 ckpt 5.1'
!***
               CALL PAH_SIZE_ATOMS(APAH,NC,NHTOT,HTOC)
               ND=DFRAC*NHTOT
               NH=NHTOT-ND
               NC=NC+418.D21*(AESVS**3-A1**3)
               NH=NC+0.25*418.D21*(AESVS**3-A1**3)
            ENDIF
         ENDIF
         NATOM=NC+NH+ND

      ELSEIF(ICASE.GE.29.AND.ICASE.LT.ICGREMT2)THEN

! DH21Ad: IDX=ICASE-29
!            = 6*JPO + JFE
!          JPO=0-9 for Poro = 0.0, 0.1, ..., 0.9
!          JFE=0-5 for f_Fe=  0.0, 0.1, ..., 0.5
! PORO  = porosity defined above
! AESVS = radius of equal-solid-volume-sphere defined above
! FVFE  = fraction of solid volume that is Fe metal defined above

         FVFE=0.12*FFE/(1.-0.22*FFE)

         FVSIL=1.-FVFE

! Astrodust density 3.41 g cm-3 as estimated by , nominal composition 
!   Mg_1.3 (Fe,Ni)_0.3 Si O_3.6
!   MW = 1.3*24+0.3*56.1+28+3.6*16=133.63
!   AW = MW/(1.3+0.3+1+3.6)=MW/6.20=21.55
!   n_atom = 9.47e22 cm-3

! Fe, density 7.87, natom=8.48e22

! Natom = V*(FVSIL*9.47e22 + FVFE*8.48e22)

         NATOM=(4.*PI/3.)*AESVS**3*(FVSIL*9.47E22+FVFE*8.48E22)

!--------------------------------------------------------------------
! following is for amorphous carbon case whenever it is implemented...
! we assume amorphous carbon density 1.85 gm cm-3
! we assume composite grain structure: up to radi us A1 is PAH-like,
! with carbon density 1.10e23 cm-3
! from A1 to A2 is graphite-like with carbon density 1.10e23 cm-3
!                            and H:C = 0.25
! beyond A2 is AMC-like with carbon density 0.925e23 cm-3
!                            and H:C = 0.25
! current settings: A1 =  50 A
!                   A2 = 100 A

!         A1=50.D-8
!         A2=100.D-8
!         IF(A.LE.A1)THEN
!            APAH=A
!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 6'
!***
!            CALL PAH_SIZE_ATOMS(APAH,NC,NHTOT,HTOC)
!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 7'
!***
!            ND=DFRAC*NHTOT
!            NH=NHTOT-ND
!         ELSE
!            APAH=A1
!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 8'
!***
!            CALL PAH_SIZE_ATOMS(APAH,NC,NHTOT,HTOC)
!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 9'
!***
!            ND=DFRAC*NHTOT
!            NH=NHTOT-ND
!            IF(A.LE.A2)THEN
!               NC=NC+460.D21*(A**3-A1**3)
!               NH=NH+0.25D0*460.D21*(A**3-A1**3)
!            ELSE
!               NC=NC+460.D21*(A2**3-A1**3)+387.D21*(A**3-A2**3)
!               NH=NH+0.25D0*(460.D21*(A2**3-A1**3)+
!     &                       387.D21*(A**3-A2**3))
!            ENDIF
!         ENDIF
!-----------------------------------------------------------------------

      ELSE
         WRITE(0,*)'vsg_td_emission_v11 Fatal error:'
         WRITE(0,*)' unknown ICASE=',ICASE
         STOP
      ENDIF

!*** diagnostic
!      write(0,*)'vsg_td_emission_v11 ckpt 10'
!***

! Have now defined NC,NH,ND for PAH, gra, or amc
!                  NATOM    for sil or Fem

! for radius AGR:
! call SIL_GRA_PAH_TEMP to obtain
!    CABS[1-NISRF] = C_abs (cm2) at wavelengths ISRF_WL[1-NISRF]
!    CSCA[1-NISRF] = C_sca (cm2) at wavelengths ISRF_WL[1-NISRF]
!    TEQ           = T(K) at which <heating> = cooling
!    EEQ           = E(erg) at temperature TEQ

! NB: SIL_GRA_PAH_TEMP also calculates for ICASE=9 (Fe metal)

!*** diagnostic
!      write(0,*)'vsg_td_emission_v11 ckpt 11'
!      write(0,*)'about to call sil_gra_pah_temp'
!      write(0,fmt='(a,1pe10.3)')' a=',a
!      write(0,*)' icase=',icase
!      write(0,fmt='(a,1pe10.3)')' natom=',natom
!      write(0,fmt='(a,1pe10.3)')' isrf(300)=',isrf(300)
!***

! Find TEQ = temperature if heating is continuous and steady

      CALL SIL_GRA_PAH_TEMP(ICASE,BOVERA,AGR,NC,NH,ND,NATOM,
     &                      DFRAC,NISRF,ISRF_WL,
     &                      ISRF,DLGLAMBDA,TEQ,
     &                      EEQ,RATE_HEATING,CABS,CSCA,SCRW,DEDT)

!*** diagnostic
      write(0,'(a)')'vsg_td_emission_v11 ckpt 12'
      write(0,'(a)')'    returned from sil_gra_pah_temp'
      write(0,fmt='(a,1pe10.3,a,1pe10.3)')'    with agr=',agr,
     &   ' TEQ=',TEQ
      write(0,fmt='(a,i3,a,1pe10.3)')'    icase=',icase,' natom=',natom
      write(0,fmt='(a,1pe10.3)')'    rate_heating=',rate_heating
!***

! if EEQ < EEQSS then seek probability distribution function
! if EEQ > EEQSS then assume steady state (1 bin only)
! try setting the cutoff at 150 eV

      EEQSS=150.0*1.60218e-12

! 17.05.15 (BTD) reduce EEQSS when TEQ is low to try to avoid
!                numerical failure

      IF(TEQ.LT.4.7D0)EEQSS=135.0*1.60218E-12

! 16.10.03 (BTD) reduce EEQSS when TEQ is low to try to avoid
!                numerical failure

      IF(TEQ.LT.4.D0)EEQSS=120.0*1.60218E-12

!*** diagnostic
!      IF(EEQ.GE.EEQSS)THEN
!         WRITE(0,*)'go with single temp *************'
!      ENDIF
!***

      IF(EEQ.LT.EEQSS)THEN

! if steady-state energy content EEQ < EEQSS, 
! calculate probability distribution for energy and temperature

! first need to set up grid of energy bins
! set UMAX for for first try: umax=13.6ev + 2*Eeq
! we require highest enthalpy bin to exceed 13.6eV

         UMAXMIN=13.65D0*1.602D-12

         IF(EEQ.LT.0.1*EEQSS)THEN
            UMIN=0.
            UMAX=13.6*1.602D-12+EEQ*2.0
         ELSE
            UMIN=EEQ/5.0
            UMAX=13.6*1.602D-12+EEQ*2.0
         ENDIF
         IF(UMAX.LT.UMAXMIN)UMAX=UMAXMIN

!*** diagnostic
!      write(0,*)'eeq(erg)=',eeq
!      write(0,*)'first try: umin,umax=',umin,umax
!***

!------------------------- carbonaceous grains ------------------------

         IF(ICASE.EQ.2.OR.     ! graphite-PAH+ as per AD12
     &      ICASE.EQ.3.OR.     ! graphite-PAH0 as per AD12
     &      ICASE.EQ.10.OR.    ! graphite-PAH+ as per HD14 a_t=0.0020um
     &      ICASE.EQ.11.OR.    ! graphite-PAH0 as per HD14 a_t=0.0020um
     &      ICASE.EQ.12.OR.    ! graphite-PAH+ as per HD14 a_t=0.0050um
     &      ICASE.EQ.13.OR.    ! graphite-PAH+ as per HD14 a_t=0.0050um
     &      ICASE.EQ.14.OR.    ! D16 graphite-PAH+ a_t=0.0020um
     &      ICASE.EQ.15.OR.    ! D16 graphite-PAH0 a_t=0.0020um
     &      ICASE.EQ.16.OR.    ! D16 graphite-PAH+ a_t=0.0050um
     &      ICASE.EQ.17.OR.    ! D16 graphite-PAH+ a_t=0.0050um
     &      ICASE.EQ.18.OR.    ! D16 graphite-PAH0 a_t=0.0020um
     &      ICASE.EQ.19.OR.    ! D16 graphite-PAH0 a_t=0.0020um
     &      ICASE.EQ.20.OR.    ! D16 graphite-PAH+ a_t=0.0050um
     &      ICASE.EQ.21.OR.    ! D16 graphite-PAH0 a_t=0.0050um
     &      ICASE.EQ.22.OR.    ! D16 graphite-PAH+ a_t=0.0020um
     &      ICASE.EQ.23.OR.    ! D16 graphite-PAH0 a_t=0.0020um
     &      ICASE.EQ.24.OR.    ! D16 graphite-PAH+ a_t=0.0050um
     &      ICASE.EQ.25.OR.    ! D16 graphite-PAH0 a_t=0.0050um
     &      ICASE.EQ.26.OR.    ! D16 graphite random c axis
     &      ICASE.EQ.27.OR.    ! D16 graphite c || a
     &      ICASE.EQ.28.OR.    ! D16 graphite, EMT  (MG,matric=E||c)
     &      ICASE.EQ.ICGREMT2)THEN  ! D16 graphite, EMT2 (MG,matrix=Eperpc)  

!*** diagnostic
!          write(0,*)'vsg_td_emission_v11 ckpt 13'
!          write(0,*)' ICASE=',ICASE
!          write(0,*)' call gra_pah_density_states ...'
!***
! Note: we are assuming that we will only be considering temperature
! fluctuations for grains small enough that NC,NH,ND < 1e9 and we therefore
! do not need to worry about overflow of integer NC_INT,NH_INT,ND_INT

            NC_INT=NINT(NC)
            NH_INT=NINT(NH)
            ND_INT=NINT(ND)

!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 14'
!            write(0,*)' about to call gra_pah_density_states with'
!            write(0,*)' NDIMMODES=',NDIMMODES
!***
            CALL GRA_PAH_DENSITY_STATES(METHOD,AGR,NC_INT,NH_INT,
     &                                  ND_INT,UMAX,NMODES,EMODES,
     &                                  JEMAX0,E,DE,LNGSUM,EBS,LNG,
     &                                  NDIME,NDIMMODES)

!              for AESVS < 2.5e-7cm: computes actual NMODES and EMODES(1-NMODES)
!                    > 2.5e-7cm; sets NMODES=50000
!                                sets EMODES(1-100) only
 
!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 15'
!            write(0,*)' returned from gra_pah_density_states with'
!            write(0,*)' de=',de
!***

!---------------------- silicate or Astrodust grains -----------------------

         ELSEIF(ICASE.EQ.1.OR.   ! astrosilicate   
     &          ICASE.EQ.4.OR.   ! HD2014 silicate, no Fe inclusion
     &          ICASE.EQ.5.OR.   ! HD2014 silicate, f_vol=0.02 Fe inclusions
     &          ICASE.EQ.6.OR.   ! HD2014 silicate, f_vol=0.05 Fe inclusions
     &          ICASE.EQ.7.OR.   ! HD2014 silicate, f_vol=0.10 Fe inclusions
     &          ICASE.EQ.8.OR.   ! HD2014 silicate, f_vol=0.20 Fe inclusions
     &          (ICASE.GT.28.AND.ICASE.LT.ICGREMT2))THEN   ! Astrodust

! NATOM has already been calculated

            NATOM_INT=NINT(NATOM)

!*** diagnostic
!            write(0,fmt='(a,i9)')
!     &         'vsg_td_emission_v11 ckpt 16, NATOM_INT=',natom_int
!***

            CALL SIL_DENSITY_STATES(METHOD,AESVS,NATOM_INT,UMAX,
     &                              NMODES,EMODES,JEMAX0,E,DE,
     &                              LNGSUM,EBS,LNG,NDIME,NDIMMODES)

!              for AESVS < 2.5e-7cm: computes actual NMODES and EMODES(1-NMODES)
!                    > 2.5e-7cm: sets NMODES=50000
!                                sets EMODES(1-100) only
!              when METHOD='dbdis', do not calculate LNG

!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 17, de=',de
!***

            IF(NATOM_INT.LT.10000)THEN
               WRITE(CDESCR(4),FMT='(I5,A)')
     &               NATOM_INT,' = N_atom for silicate'
            ELSE
               WRITE(CDESCR(4),FMT='(1PE10.3,A)')
     &               NATOM,' = N_atom for silicate'
            ENDIF
         ELSEIF(ICASE.EQ.9)THEN

! ICASE=9 : metallic Fe

            NATOM_INT=NINT(NATOM)

!*** diagnostic
!            write(0,'(a,i6)')'vsg_td_emission_v11 ckpt 18, natom_int=',
!     &         natom_int
!***
            CALL FE_DENSITY_STATES(METHOD,AESVS,NATOM_INT,UMAX,
     &                             NMODES,EMODES,JEMAX0,E,DE,
     &                             LNGSUM,EBS,LNG,NDIME,NDIMMODES)

!              for A < 2.5e-7cm: computes actual NMODES and EMODES(1-NMODES)
!                    > 2.5e-7cm; sets NMODES=50000
!                                sets EMODES(1-100) only
!              when METHOD='dbdis', do not calculate LNG 

!*** diagnostic
!            write(0,'(a,i6,a,1pe10.3)')
!     &         'vsg_td_emission_v11 ckpt 19, natom_int=',
!     &         natom_int,'de=',de
!***
         ELSE
            WRITE(0,*)'vsg_td_emission_v11 ckpt 20'
            WRITE(0,*)'Fatal error: unrecognized ICASE=',ICASE
            STOP
         ENDIF

!*** diagnostic
!         write(0,*)'vsg_td_emission_v11 ckpt 21, icase=',icase
!***

         WRITE(CDESCR(5),FMT='(A)')' '

! start iteration to find solution PSTATE(1-NSTAT)

 0600    CONTINUE   ! re-enter here if have adjusted UMIN,UMAX
                    ! recalculate energy bins U(1-NSTATE) for new UMIN,UMAX
                    ! recalculate PSTATE

! -------------------------------------------------------------------  
! construct enthalpy bins; 
! input: umax (the highest bin of the enthalpy space); 
!        nstate -- the number of states;
! output: u(1-nstate) (enthalpy, erg); 
!        ua(1-nstate) (lower boundary, erg) 
!        ub(1-nstate) (upper boundary, erg)
!        lngu(1-nstate) (density of states)
!        EEQ = equilibrium enthalpy (erg)

         DO I=1,NSTATE
            U(I)=0.D0
            UA(I)=0.D0
            UB(I)=0.D0
            T(I)=0.D0
            DELTA_T(I)=0.D0
         ENDDO
	
!*** diagnostic
!         write(0,fmt='(a,i3,a,1pe10.3)')
!     &      'vsg_td_emission_v11 ckpt 22, icase=',icase,' natom=',natom
!***

!*** diagnostic
         WRITE(0,7601)UMIN/HC,UMAX/HC
!***

         IF(ICASE.EQ.2.OR.   ! gra_PAHi_D2003: gra-PAH+
     &      ICASE.EQ.3.OR.   ! gra_PAHn_D2003: gra-PAH0
     &      ICASE.EQ.10.OR.  ! gra_PAHi__2014: gra-PAH+ a_t=0.0020um
     &      ICASE.EQ.11.OR.  ! gra_PAHn__2014: gra-PAH0 a_t=0.0020um
     &      ICASE.EQ.12.OR.  ! gra_PAHib_2014: gra-PAH+ a_t=0.0050um
     &      ICASE.EQ.13.OR.  ! gra_PAHnb_2014: gra-PAH0 a_t=0.0050um
     &      ICASE.EQ.14.OR.  ! gra_D16rc_PAHi: D16gra,rand c, PAH+ a_t=0.0020um
     &      ICASE.EQ.15.OR.  ! gra_D16rc_PAHn: D16gra,rand c, PAH0 a_t=0.0020um
     &      ICASE.EQ.16.OR.  ! gra_D16rcPAHib: D16gra,rand c, PAH+ a_t=0.0050um
     &      ICASE.EQ.17.OR.  ! gra_D16rcPAHnb: D16gra,rand c, PAH0 a_t=0.0050um
     &      ICASE.EQ.18.OR.  ! gra_D16ca_PAHi: D16gra,c || a, PAH+ a_t=0.0020um
     &      ICASE.EQ.19.OR.  ! gra_D16ca_PAHn: D16gra,c || a, PAH0 a_t=0.0020um
     &      ICASE.EQ.20.OR.  ! gra_D16caPAHib: D16gra,c || a, PAH+ a_t=0.0050um
     &      ICASE.EQ.21.OR.  ! gra_D16caPAHnb: D16gra,c || a, PAH0 a_t=0.0050um
     &      ICASE.EQ.22.OR.  ! gra_D16emtPAHi: D16gra, EMT, PAH+ a_t=0.0020um
     &      ICASE.EQ.23.OR.  ! gra_D16emtPAHn: D16gra, EMT, PAH0 a_t=0.0020um
     &      ICASE.EQ.24.OR.  ! graD16emtPAHib: D16gra, EMT, PAH+ a_t=0.0050um
     &      ICASE.EQ.25.OR.  ! graD16emtPAHnb: D16gra, EMT, PAH0 a_t=0.0050um
     &      ICASE.EQ.26.OR.  ! graphite_D16rc: D16 graphite, random c axis
     &      ICASE.EQ.27.OR.  ! graphite_D16ca: D16 graphite, c || a
     &      ICASE.EQ.28.OR.  ! graphiteD16emt: D16 graphite, EMT
     &      ICASE.EQ.ICGREMT2.OR. ! D16 graphite, EMT (MG2; matrix=Eperpc)
     &      ICASE.EQ.90.OR.ICASE.EQ.91)THEN ! D16 EMT2
!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 23'
!            write(0,*)'call gra_pah_enthalpy_state with METHOD=',METHOD
!            write(0,*)'                                 DE=',DE
!***
            CALL GRA_PAH_ENTHALPY_STATE(NC_INT,NH_INT,ND_INT,UMIN,UMAX,
     &                                  NSTATE,U,UA,UB,METHOD,T,DELTA_T,
     &                                  NSET,NMODES,EMODES,DE)

!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 24'
!            write(0,*)'returned from gra_pah_enthalpy_state...'
!            write(0,*)'call gra_pah_degeneracy...'
!***
            CALL GRA_PAH_DEGENERACY(METHOD,AESVS,NC_INT,NH_INT,ND_INT,
     &                              UMAX,NSTATE,U,UA,UB,LNGU,NSET,
     &                              JEMAX0,E,DE,LNGSUM,EBS)

!              if METHOD='stati', calculate LNGU(1-NSTATE)
!              if METHOD='dbdis', do not calculate LNGU [set = 0] 

!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 25'
!            write(0,*)'returned from gra_pah_degeneracy...'
!***

         ELSEIF(ICASE.EQ.1.OR.
     &          ICASE.EQ.4.OR.
     &          ICASE.EQ.5.OR.
     &          ICASE.EQ.6.OR.
     &          ICASE.EQ.7.OR.
     &          ICASE.EQ.8.OR.
     &          (ICASE.GE.29.AND.ICASE.LT.ICGREMT2))THEN   ! Astrodust

!*** diagnostic
            write(0,fmt='(a,/,a,1pe10.3,a)')
     &         'vsg_td_emission_v11 ckpt 26:',
     &         '   call sil_enthalpy_state with de=',de,
     &         ' cm-1 for Beyer-Swinehart alg.'
!***

! for the moment, assume density of states is that of silicate
! even when Fe inclusions are present

            CALL SIL_ENTHALPY_STATE(AESVS,UMIN,UMAX,NSTATE,U,UA,UB,
     &                              METHOD,T,DELTA_T,NSET,NMODES,
     &                              EMODES,DE)

!              calculate U(1-NSTATE) between UMIN and UMAX
!                        UA(1-NSTATE),UB(1-NSTATE)
!              if Umin=0, first 3 bins correspond to modes EMODES(1-3)
!                         bins 4-NSET are chosen to sample low energy
!                         mode spectrum
!                         bins NSET+1-NSTATE are sensibly distributed
!                         (NSET .le. 11)

!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 27'
!***
            CALL SIL_DEGENERACY(METHOD,AESVS,UMAX,NSTATE,
     &                          U,UA,UB,LNGU,NSET,JEMAX0,E,DE,
     &                          LNGSUM,EBS)

!              if METHOD='stati', calculate LNGU(1-NSTATE)
!              if METHOD='dbdis', do not calculate LNGU [set = 0] 

!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 28'
!***

         ELSEIF(ICASE.EQ.9)THEN

!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 29: icase=',icase
!            write(0,fmt='(a,1pe10.3)')
!     &         'about to call fe_enthalpy_state with de=',de
!***
            CALL FE_ENTHALPY_STATE(AESVS,UMIN,UMAX,NSTATE,U,UA,UB,
     &                             METHOD,T,DELTA_T,NSET,NMODES,
     &                             EMODES,DE)

!              calculate U(1-NSTATE) between UMIN and UMAX
!                        UA(1-NSTATE),UB(1-NSTATE)
!              if Umin=0, first 3 bins correspond to modes EMODES(1-3)
!                         bins 4-NSET are chosen to sample low energy
!                         mode spectrum
!                         bins NSET+1-NSTATE are sensibly distributed
!                         (NSET .le. 11)

!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 30:'
!            write(0,fmt='(a,i5)')'   nstate=',nstate
!            write(0,fmt='(a,1pe10.3)')'   u(nstate)=',u(nstate)
!            write(0,fmt='(a,1pe10.3)')'   T(nstate)=',T(nstate)
!            write(0,*)'about to call fe_degeneracy'
!***
            CALL FE_DEGENERACY(METHOD,AESVS,UMAX,NSTATE,
     &                          U,UA,UB,LNGU,NSET,JEMAX0,E,DE,
     &                          LNGSUM,EBS)

!              if METHOD='stati', calculate LNGU(1-NSTATE)
!              if METHOD='dbdis', do not calculate LNGU [set = 0] 

!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 31'
!***
         ELSE
            WRITE(0,*)'vsg_td_emission_v11 ckpt 32:'
            WRITE(0,*)' unrecognized icase=',icase
            STOP
         ENDIF

!*** diagnostic
!         write(0,*)'vsg_td_emission_v11 ckpt 33: icase=',icase
!***
! initialize 
!    Amatrix(1:nstate,1:nstate),
!    Pstate(nstate),
!    Amatrix1(1:nstate1,1:nstate1),
!    Bmatrix(1:nstate1),
!    rhs(1:nstate1), 
!    Pstate1(1:nstate1), 
! and ija(1:nmax),sa(1:nmax). 
! modified on 99.11.28.

         DO I=1,NSTATE
            PSTATE(I)=0.D0
	    PCUMULT(I)=0.D0
            COOLING_TIME(I)=0.D0
            RAD_TIME(I)=0.D0
            DO F=1,NSTATE
               AMATRIX(F,I)=0.D0
               BGD(F,I)=0.D0
            ENDDO
         ENDDO

!*** diagnostic
!         write(0,*)'vsg_td_emission_v11 ckpt 33.5'
!         write(0,*)'           NSTATE=',NSTATE,' NSTATE1=',NSTATE1
!*** 
         DO I=1,NSTATE1
            PSTATE1(I)=0.D0
            RHS(I)=0.D0
!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 33.6, I=',I
!***
            DO F=1,NSTATE1
!*** diagnostic
!               write(0,*)'vsg_td_emission_v11 ckpt 33.7, I,F=',I,F
!***
               AMATRIX1(F,I)=0.D0
            ENDDO
         ENDDO

! construct the transition matrix;

!*** diagnostic
!      write(0,*)'vsg_td_emission_v11 ckpt 34 :'
!      write(0,*)' call transition_matrix with'
!      write(0,*)' method,icase,ishap=',method,icase,ishap
!      write(0,*)' a=',a
!      write(0,*)' dfrac=',dfrac
!***

         CALL TRANSITION_MATRIX(METHOD,ICASE,BOVERA,AGR,
     &                          DFRAC,NISRF,ISRF_WL,ISRF,DLGLAMBDA,
     &                          CABS,CSCA,NSTATE,U,UA,UB,LNGU,T,
     &                          AMATRIX,COOLING_TIME,RAD_TIME,
     &                          PSTATE,BGD)

!*** diagnostic
!      write(0,*)'vsg_td_emission_v11 ckpt 35 :'
!      write(0,*)' returned from transition_matrix'
!***

! -------------------------------------------------------------------

         WRITE(0,fmt='(a,a)')'vsg_td_emission_v11 ckpt 36: ',
     &                       'Solve the N*N matrix (BiCG) !'
         WRITE(0,fmt='(a,i4)')'   NSTATE1=',NSTATE1

! solve the linear equations: Amatrix[,]*Pstate=rhs[nstate]; 	
! we adopt BiCG method; 
! GAUSS, GMRES, BiCGstab(l<=2), BiCGstab (l>2) have been tried;
! BiCG is good enough; Gauss method not good;
! rhs: righ-hand-size vector 
! (Amatrix1[,]*Pstate1[nstate1]=rhs[nstate1]);

! calculate rate of transitions out of bin with U=8.D-12 erg = 5 eV
! and compare to total photoabsorption rate.
! if photoabsorption rate is sufficiently low, then obtain T distribution
! function using the low-U approximation

         F=1
         DO I=1,NSTATE
            IF(U(I).LT.8.D-12)F=I
         ENDDO

! AMATRIX(1,1) = total rate out of level 1
!              = photoabsorption rate
! AMATRIX(F,F) = total rate out of level F 
!                (photoabsorption plus spontaneous emission)

! TERM = (photoabsorption rate)/(total rate out of 5 eV bin)

         TERM=AMATRIX(1,1)/AMATRIX(F,F)

!*** diagnostic
         WRITE(0,fmt='(a,/,1pe12.3,a,/,1pe12.3,a)')
     &      'vsg_td_emission_v11 ckpt 37:',
     &      -AMATRIX(1,1),' =-AMATRIX(1,1)= total rate out of level 1',
     &      TERM,' =TERM= (photoabs rate)/(total rate out of 5eV bin)'
!***

! 2017.06.07 (BTD) trial value: 1.d-8

         IF(TERM.LT.1.D-8)THEN

! ------------------- low U approximation for dP/dT ------------------

!*** diagnostic
            WRITE(0,fmt='(a)')
     &         'vsg_td_emission_v11 ckpt 38: low U approx'
!***
            PSTATE(1)=1.D0
            PSTATE(NSTATE)=-AMATRIX(NSTATE,1)/AMATRIX(NSTATE,NSTATE)
            DO F=NSTATE-1,2,-1
               TERM=AMATRIX(F,1)
               DO I=F+1,NSTATE
                  TERM=TERM+PSTATE(I)*AMATRIX(F,I)
               ENDDO
               PSTATE(F)=-TERM/AMATRIX(F,F)
            ENDDO

! normalize PSTATE

            TERM=PSTATE(1)
            DO I=2,NSTATE
               TERM=TERM+PSTATE(I)
            ENDDO
            DO I=1,NSTATE
               PSTATE(I)=PSTATE(I)/TERM
            ENDDO
            PMAX=PSTATE(1)

            E_ABS=0.D0
            E_EM=0.D0
            DO I=1,NSTATE1-1
               DO F=I+1,NSTATE1
                  E_ABS=E_ABS+PSTATE(I)*AMATRIX(F,I)*(U(F)-U(I))
               ENDDO
            ENDDO
            DO I=2,NSTATE1
               DO F=1,I-1
                  E_EM=E_EM+PSTATE(I)*AMATRIX(F,I)*(U(I)-U(F))
               ENDDO
            ENDDO
            DO I=1,NSTATE1
               PSTATE1(I)=PSTATE(I)
            ENDDO

! E_EM,E_ABS are in erg s-1 (per particle)

            write(0,9876)1,PSTATE1(1)/PMAX,U(1)/hc,T(1)
            write(0,9876)2,PSTATE1(2)/PMAX,U(2)/hc,T(2)
            write(0,9876)3,PSTATE1(3)/PMAX,U(3)/hc,T(3)
            write(0,9876)4,PSTATE1(4)/PMAX,U(4)/hc,T(4)
            write(0,9876)NSTATE1-3,PSTATE1(NSTATE1-3)/PMAX,U(NSTATE1-3)/
     &                   hc,T(NSTATE1-3)
            write(0,9876)NSTATE1-2,PSTATE1(NSTATE1-2)/PMAX,U(NSTATE1-2)/
     &                   hc,T(NSTATE1-2)
            write(0,9876)NSTATE1-1,PSTATE1(NSTATE1-1)/PMAX,U(NSTATE1-1)/
     &                   hc,T(NSTATE1-1)
            write(0,9876)NSTATE1,PSTATE1(NSTATE1)/PMAX,U(NSTATE1)/
     &                   hc,T(NSTATE1)
            WRITE(0,7901)E_EM,E_ABS,(E_EM/E_ABS-1.D0)               

         ELSE

!*** diagnostic
            WRITE(0,fmt='(a)')
     &         'vsg_td_emission_v11 ckpt 38.5: full solution'
!***

!-----------------------------------------------------------------------
!
!                   N              N
! 0 = (d/dt) P_i = sum A_ij P_j - sum A_ki P_i
!                 j.neq.i        k.neq.i
!
!                  N
! Define A_ii = - sum A_ki
!                k.neq.i
!
! Then
!      N
! 0 = sum A_ij P_j   for i = 1-N
!     j=1
!
! with normalization
!      N
! 1 = sum P_j
!     j=1
!                              N-1
! Thus, we can write P_N = 1 - sum P_j
!                              j=1
!
! and obtain
!     N-1                        N-1
! 0 = sum A_ij P_j + A_iN - A_iN sum P_j
!     j=1                        j=1
!
! Define new matrix B_ij  for i=1,...,N-1, j=1,...,N-1:
!
! B_ij = A_ij - A_iN     [B_ij = matrix AMATRIX1(i,j) below]
!
! Then we have N-1 equations for the unknown P_j, j=1,...,N-1
! N-1
! sum B_ij P_j = -A_iN   [-A_iN = RHS(i) below]
! j=1
!
!-----------------------------------------------------------------------

            DO F=1,NSTATE1
               RHS(F)=-AMATRIX(F,NSTATE)
               PSTATE1(F)=0.
               DO I=1,NSTATE1
                  AMATRIX1(F,I)=AMATRIX(F,I)-AMATRIX(F,NSTATE)
               ENDDO
            ENDDO

! initialize ija[1:nmax] and sa[1:nmax];

            DO I=1,NMAX
               IJA(I)=0
               SA(I)=0.0
            ENDDO

            THRESH=1.D-90
!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 39'
!*** 
            CALL DSPRSIN(AMATRIX1,NSTATE1,THRESH,NMAX,SA,IJA)
!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 40'
!***
            ITOL=3

! 2006.11.04 try smaller TOL
!      TOL=1.0D-14

            TOL=1.0D-15

            IT_END=0

            ITMAX=1000
 1100       ITMAX=ITMAX+100
!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 41'
!***
            CALL LINBCG(NMAX,SA,IJA,NSTATE1,RHS,PSTATE1,
     &                  ITOL,TOL,ITMAX,ITER,ERR)

!*** diagnostic
!            write(0,*)'vsg_td_emission_v11 ckpt 42'
!***
            IT_END=IT_END+ITER
            IF(ERR.GT.TOL)GOTO 1100
            WRITE(0,7700)IT_END,ERR,TOL

!-----------------------------------------------------------------------
! 2011.05.19 (BTD) check that solution has actually been found
!            if error is unacceptable, then return to LINGCG for
!            another round...

            WRITE(0,FMT='(A)')
     &        'vsg_td_emission_v11 ckpt 43: check solution'
            REFINE=.FALSE.
            DO F=1,NSTATE1
               TERM=-RHS(F)
               DO I=1,NSTATE1
                  TERM=TERM+AMATRIX1(F,I)*PSTATE1(I)
               ENDDO
               TERM=ABS(TERM/AMATRIX1(1,1))
               IF(PSTATE1(F).GT.1.D-10.AND.TERM.GT.1.D-6)THEN
                  WRITE(0,7900)F,PSTATE1(F),TERM
                  REFINE=.TRUE.
               ENDIF
            ENDDO
            E_ABS=0.D0
            E_EM=0.D0
            DO I=1,NSTATE1-1
               DO F=I+1,NSTATE1
                  E_ABS=E_ABS+PSTATE1(I)*AMATRIX(F,I)*(U(F)-U(I))
               ENDDO
            ENDDO
            DO I=2,NSTATE1
               DO F=1,I-1
                  E_EM=E_EM+PSTATE1(I)*AMATRIX(F,I)*(U(I)-U(F))
               ENDDO
            ENDDO

c E_EM,E_ABS are in erg s-1 (per particle)

            WRITE(0,7901)E_EM,E_ABS,(E_EM/E_ABS-1.D0)
            IF(REFINE)GOTO 1100

!-----------------------------------------------------------------------
            SUM_PSTATE=0.

!*** diagnostic
!         write(0,*)'vsg_td_emission_v11 ckpt 44: nstate1=',nstate1
!***
            DO I=1,NSTATE1

! 080308: add check for possible numerical failure in LINBCG:

               TERM=ABS(PSTATE1(I))
               IF(.NOT.(ABS(PSTATE1(I))-TERM).EQ.0.)THEN
                  WRITE(0,*)'I=',I,' PSTATE1(I)=',PSTATE1(I)
                  WRITE(0,*)'Fatal error: STOP'
                  STOP
               ENDIF

! it is possible that some level pops may be small negative numbers
! if so, set them to zero

               IF(PSTATE1(I).LT.0.D0)PSTATE1(I)=0.D0
               SUM_PSTATE=SUM_PSTATE+PSTATE1(I)
            ENDDO	
            WRITE(0,7701)SUM_PSTATE

! ------------------------------------------------------------------------

! fine tuning the temperature range [umin, umax] 
! for iteration; modified on Jan.27, 2000: 
! umin fixed at 0.;

            PMAX=0.
            DO I=1,NSTATE1
               IF(PSTATE1(I).GE.PMAX)THEN
                  PMAX=PSTATE1(I)
                  JPMAX=I
               ENDIF
            ENDDO

!*** diagnostic
            write(0,7610)PMAX,JPMAX,U(JPMAX)/hc,T(JPMAX)
!***

!-----------------------------------------------------------------------
!                         new module for adjusting Umax
! 2006.11.04 reduce PMIN_LO and PMIN_UP
!      PMIN_LO=1.D-11
!      PMIN_UP=1.D-11
            PMIN_LO=1.D-13
            PMIN_UP=1.D-13
            REFINE=.FALSE.

!*** diagnostic
            write(0,9876)1,PSTATE1(1)/PMAX,U(1)/hc,T(1)
            write(0,9876)2,PSTATE1(2)/PMAX,U(2)/hc,T(2)
            write(0,9876)3,PSTATE1(3)/PMAX,U(3)/hc,T(3)
            write(0,9876)4,PSTATE1(4)/PMAX,U(4)/hc,T(4)
            write(0,9876)NSTATE1-3,PSTATE1(NSTATE1-3)/PMAX,U(NSTATE1-3)/
     &                   hc,T(NSTATE1-3)
            write(0,9876)NSTATE1-2,PSTATE1(NSTATE1-2)/PMAX,U(NSTATE1-2)/
     &                   hc,T(NSTATE1-2)
            write(0,9876)NSTATE1-1,PSTATE1(NSTATE1-1)/PMAX,U(NSTATE1-1)/
     &                   hc,T(NSTATE1-1)
            write(0,9876)NSTATE1,PSTATE1(NSTATE1)/PMAX,U(NSTATE1)/
     &                   hc,T(NSTATE1)
!*** 

! ---------- adjust UMAX if necessary -----------------------------------

            IF(PSTATE1(NSTATE1)/PMAX.LE.PMIN_UP.AND.
     &         UMAX.GT.UMAXMIN)THEN

! Need to decrease UMAX
! find point JCUT where we could cut off distribution
! start from highest bin and move down to first bin with P/Pmax > Pmin

               write(0,7620)umax/hc

               UMAXHI=UMAX
               DO I=NSTATE1,1,-1
                  JCUT=I
                  IF(UB(I).LT.UMAXMIN)GOTO 2200
                  IF(PSTATE1(I)/PMAX.GT.PMIN_UP)GOTO 2200
               ENDDO

! check that there are no bins below this with P.le.0

 2200          DO I=JPMAX,JCUT
                  IF(PSTATE1(I).LE.0.)THEN
                     JCUT=I

!*** diagnostic
                     write(0,*)'*** reduced jcut to jcut=',jcut
!***

                     GOTO 2400
                  ENDIF
               ENDDO

! if UMAX is less than 1.02*UMAXLO, stop iterating
! 2006.11.04 make termination condition more stringent :1.05 -> 1.02
! if UMAX is less than 1.05*UB(JCUT), stop iterating
! if UMAX is less than 1.01*UB(JCUT), stop iterating
! otherwise, reduce UMAX and try again

 2400          IF(UMAX.GT.1.02*UMAXLO.AND.UMAX.GT.1.01*UB(JCUT))THEN
                  UMAX=0.8*UB(JCUT)+0.2*UMAX
                  IF(UMAX.LT.UMAXLO)UMAX=1.01*UMAXLO
                  IF(UMAX.LT.UMAXMIN)UMAX=UMAXMIN
                  REFINE=.TRUE.

!*** diagnostic
                  write(0,7650)JCUT,UMAX/HC
!***
               ENDIF

            ELSEIF(PSTATE1(NSTATE1)/PMAX.GT.PMIN_UP)THEN

! need to increase UMAX:

               UMAXLO=UMAX
!---------new
               FAC=1.2
               IF(PSTATE1(NSTATE1)/PMAX.LT.1.E7*PMIN_UP)FAC=1.1
               IF(PSTATE1(NSTATE1)/PMAX.LT.1.E2*PMIN_UP)FAC=1.05
               IF(FAC*UMAX.LT.UMAXHI)THEN
                  REFINE=.TRUE.
                  UMAX=FAC*UMAX
                  WRITE(0,7670)UMAXHI/HC,FAC,UMAX/HC

!---------
!               IF(1.2*UMAX.LT.UMAXHI)THEN
!
! 2020.05.08 (BTD) change to avoid too large jump in Umax
!
! 2020.05.14 (BTD) experiment
!                  IF(PSTATE1(NSTATE1)/PMAX.LT.1.E2*PMIN_UP)THEN
!                  IF(PSTATE1(NSTATE1)/PMAX.LT.1.E7*PMIN_UP)THEN
!----------------------------
!                     UMAX=1.1*UMAX
!                     write(0,7671)umaxhi/hc,umax/hc
!                  ELSE
!                     UMAX=1.2*UMAX
!                     write(0,7670)umaxhi/hc,umax/hc
!                  ENDIF
!                 REFINE=.TRUE.
!---------------------------------------------------------------
               ELSE
                  IF(UMAX/UMAXHI-1.GT.0.01)THEN
                     UMAX=.5*(UMAXHI+UMAX)
                     REFINE=.TRUE.
                     write(0,7680)umaxhi/hc,umax/hc
                  ENDIF
               ENDIF

            ENDIF

! --------------- adjust UMIN if necessary ----------------------------

            IF(PSTATE1(2)/PMAX.LT.PMIN_LO)THEN

! need to increase Umin

               UMINLO=UMIN

               DO I=1,NSTATE1
                  JCUT=I
                  IF(PSTATE1(I)/PMAX.GT.PMIN_LO)GOTO 2600
               ENDDO
 2600          IF(UMIN.LT.0.95*UA(JCUT))THEN
                  UMIN=.2*UMIN+.8*UA(JCUT)
                  REFINE=.TRUE.
                  write(0,7690)JCUT,uminlo/hc,umin/hc
               ENDIF

! 2007.12.27 (BTD) change: do not allow UMIN/hc to be decreased below 1.0 cm-1

            ELSEIF(UMIN.GT.HC.AND.PSTATE1(1)/PMAX.GT.PMIN_LO)THEN

            write(0,*)'--- attempt to reduce umin ---'

! need to reduce Umin

               UMINHI=UMIN
               IF(0.8*UMIN.GT.UMINLO)THEN

! 2020.05.14 (BTD) change to avoid too large reduction in Umin

                  IF(PSTATE1(1)/PMAX.LT.1.E2*PMIN_LO)THEN
                     UMIN=MAX(HC,0.95*UMIN)
                     WRITE(0,7694)UMINLO/HC,UMIN/HC
                  ELSE
                     IF(PSTATE1(1)/PMAX.LT.1.E3*PMIN_LO)THEN
                        UMIN=MAX(HC,0.90*UMIN)
                        WRITE(0,7695)UMINLO/HC,UMIN/HC
                     ELSE
                        UMIN=MAX(HC,0.8*UMIN)
                        WRITE(0,7696)UMINLO/HC,UMIN/HC
                     ENDIF
                  ENDIF
                  REFINE=.TRUE.
               ELSE
                  IF(UMIN-UMINLO.GT.0.01*UMIN.AND.
     &               UMIN/HC.GT.20.)THEN
                     UMIN=MAX(0.5*(UMINLO+UMIN),HC)
                     REFINE=.TRUE.
                  ENDIF
               ENDIF
            ENDIF

            IF(REFINE)GOTO 0600
            WRITE(0,*)'Iteration complete:'
         ENDIF
         WRITE(0,*)'Store dP/dT and compute emitted spectrum'

!------------------------------------------------------------------------
! end iteration;

         PCUMULT(NSTATE1)=PSTATE1(NSTATE1)
         DO I=NSTATE1-1,1,-1
            PCUMULT(I)=PCUMULT(I+1)+PSTATE1(I)
         ENDDO

      ELSE

! delta-function T distribution

         NSTATE1=1
         PSTATE1(1)=1.
         PCUMULT(1)=1.
         T(1)=TEQ
         U(1)=EEQ
         UB(1)=1.0005*EEQ
         UA(1)=0.9995*EEQ
         COOLING_TIME(1)=EEQ/RATE_HEATING
         RAD_TIME(1)=COOLING_TIME(1)

!*** diagnostic
!         write(0,*)'vsg_td_emission_v11 ckpt 45'
!***
      ENDIF

!*** diagnostic
!      write(0,fmt='(a,1pe10.3)')'vsg_td_emission_v11 ckpt 46, a=',a
!***
!-----------------------------------------------------------------------
! ------------------  write output to files  ---------------------------
! unit=66 = vsg_stat_therm.pout
!      70 = vsg_stat_therm.dpdtlib.out

!      WRITE(66,6601)(CDESCR(I),I=1,9)
!      WRITE(70,7001)NSTATE1,1.E8*A,CASEDESCRIPT,RADDESCRIPT,RATE_HEATING
!      WRITE(70,7002)
!      DO I=1,NSTATE1
!         IF(PSTATE1(I).GE.0.D0)THEN
!            IF(T(I).LT.1.E4)THEN
!               WRITE(66,6660)T(I),U(I)/HC,(UB(I)-UA(I))/HC,LNGU(I),
!     &                       PSTATE1(I),PCUMULT(I),COOLING_TIME(I),
!     &                       RAD_TIME(I)
!            ELSE
!               WRITE(66,6670)T(I),U(I)/HC,(UB(I)-UA(I))/HC,LNGU(I),
!     &                       PSTATE1(I),PCUMULT(I),COOLING_TIME(I),
!     &                       RAD_TIME(I)
!            ENDIF
!         ELSE
!            IF(T(I).LT.1.E4)THEN
!               WRITE(66,6661)T(I),U(I)/HC,(UB(I)-UA(I))/HC,LNGU(I),
!     &                       PSTATE1(I),PCUMULT(I),COOLING_TIME(I),
!     &                       RAD_TIME(I)
!            ELSE
!               WRITE(66,6671)T(I),U(I)/HC,(UB(I)-UA(I))/HC,LNGU(I),
!     &                       PSTATE1(I),PCUMULT(I),COOLING_TIME(I),
!     &                       RAD_TIME(I)
!            ENDIF
!         ENDIF
!         IF(T(I).LT.1.E4)THEN
!            WRITE(70,7010)T(I),U(I)/HC,(UB(I)-UA(I))/HC,PSTATE1(I)
!         ELSE
!            WRITE(70,7011)T(I),U(I)/HC,(UB(I)-UA(I))/HC,PSTATE1(I)
!         ENDIF
!      ENDDO

      PSTATE(NSTATE)=0.D0
      DO I=1,NSTATE1
         PSTATE(I)=PSTATE1(I)
         PSTATE(NSTATE)=PSTATE(NSTATE)+AMATRIX(NSTATE,I)*PSTATE1(I)
      ENDDO
      PSTATE(NSTATE)=-PSTATE(NSTATE)/AMATRIX(NSTATE,NSTATE)

! calculate the emission spectrum for a given size:

      IF(IWRITE.EQ.1)WRITE(77,7710)(CDESCR(I),I=1,9)
      PRAD=0.
      DO I=1,NISRF
         EMISSION(I)=0.

! if we have a probability distribution, calculate emission

         IF(NSTATE1.GT.1)THEN
            DO J=2,NSTATE1
               EMSN_INTRABIN=0.D0
               EMSN_INTERBIN=0.D0

               IF(METHOD.EQ.'stati')THEN

                  IF((H*C/ISRF_WL(I)).LE.(UB(J)-UA(J)))THEN
                     EMSN_INTRABIN=(1.+ISRF(I)/(8.*PI*HCC/
     &                              ISRF_WL(I)**5))*
     &                             (8.*PI*HCC/(ISRF_WL(I)**4))*
     &                             CABS(I)*PSTATE(J)*
     &                             (1.-H*C/ISRF_WL(I)/(UB(J)-UA(J)))
                  ENDIF

                  DO K=J-1,1,-1
                     CALL GLU(ISRF_WL(I),UA(J),UB(J),UA(K),UB(K),GKJ)
                     IF((H*C/ISRF_WL(I)).LE.(UB(J)-UA(K)).AND.
     &                  (H*C/ISRF_WL(I)).GE.(UA(J)-UB(K)))THEN
                        EMSN_INTERBIN=EMSN_INTERBIN+
     &                     (1.+ISRF(I)/(8.*PI*HCC/ISRF_WL(I)**5))*
     &                     (8.*PI*HCC/(ISRF_WL(I)**4))*CABS(I)*
     &                     PSTATE(J)*DEXP(LNGU(K)-LNGU(J))*GKJ

                     ENDIF
                  ENDDO

               ELSEIF(METHOD.NE.'stati')THEN

                  IF((H*C/ISRF_WL(I)).LE.UB(J))THEN
                     DJ=(1.+ISRF(I)/(8.*PI*HCC/ISRF_WL(I)**5))*
     &                  (8.*PI*HCC/ISRF_WL(I)**4)*CABS(I)*
     &                  PSTATE(J)/(DEXP(H*C/(ISRF_WL(I)*KB*T(J)))-1.D0)
                     EMSN_INTERBIN=EMSN_INTERBIN+DJ
                  ENDIF
               ENDIF

! calculate EMISSION=nu*P_nu, where 
! P_nu*dnu = power radiated by grain into frequency interval dnu
! nu*P_nu  = power per unit ln(lambda) or ln(nu)

               EMISSION(I)=EMISSION(I)+EMSN_INTRABIN+EMSN_INTERBIN
            ENDDO

         ELSE

! calculate thermal emission at temperature TEQ

            TERM=H*C/(KB*TEQ*ISRF_WL(I))
            IF(TERM.LT.50.)THEN
               EMISSION(I)=(1.+ISRF(I)/(8.*PI*HCC/ISRF_WL(I)**5))*
     &                     (8.*PI*HCC/ISRF_WL(I)**4)*CABS(I)/
     &                     (DEXP(TERM)-1.D0)
            ELSE
               EMISSION(I)=0.
            ENDIF
         ENDIF

! if nu*P_nu < 1e-99 erg s-1, set to zero

         IF(EMISSION(I).LT.1.D-99)EMISSION(I)=0.

         PRAD=PRAD+EMISSION(I)

         IF(IWRITE.EQ.1)WRITE(77,7750)ISRF_WL(I)*1.D4,EMISSION(I),
     &                                CABS(I)

      ENDDO
      PRAD=PRAD*LOG(ISRF_WL(11)/ISRF_WL(1))/10.
      IF(NSTATE1.GT.1)THEN
         WRITE(0,7770)PRAD
      ELSE
         WRITE(0,7771)TEQ,PRAD
      ENDIF

!*** diagnostic
!      write(0,fmt='(a,i5,a,1pe10.3)')
!     &   'vsg_td_emission_v11 ckpt 47, NSTATE1=',NSTATE1,' a=',a
!***

! When distribution function has been found (as opposed to delta fn)
! check energy balance in terms of matrix transitions;
! E_abs = E_em ?
! ==============

      IF(NSTATE1.GT.1)THEN
         E_ABS=0.
         DO I=1,NSTATE1
            DO J=I+1,NSTATE
               E_ABS=E_ABS+AMATRIX(J,I)*PSTATE(I)*(U(J)-U(I))
            ENDDO
         ENDDO
	
         E_EM=0.
         DO I=2,NSTATE
            DO J=1,I-1	
               E_EM=E_EM+AMATRIX(J,I)*PSTATE(I)*(U(I)-U(J))
            ENDDO
         ENDDO

! calculate contribution of the highest state to the emission

         TERM=0.
         DO J=1,NSTATE-1
            TERM=TERM+AMATRIX(J,NSTATE)*PSTATE(NSTATE)*(U(NSTATE)-U(J))
         ENDDO
         WRITE(0,7780)AGR*1.D8,E_ABS,E_EM,TERM

! very simple (but quite accurate!) estimation of 
! absorption/emision rate; accuracies were checked 
! by Simpson integration method;

         DO I=1,NISRF
            ABSORPTION(I)=CABS(I)*ISRF(I)
         ENDDO

         CALL TRAP_INTEG(NISRF,ISRF_WL,ABSORPTION,DLGLAMBDA,E_ABS)

         DELTA_E=DABS(E_ABS-E_EM)/E_ABS

         WRITE(0,7800)E_ABS,DELTA_E

! sanity check: halt if energy conservation fails by more than 3%

         IF(DELTA_E.GT.0.03D0)THEN
            WRITE(0,*)'vsg_td_emission_v11 ckpt 99'
            WRITE(0,*)'********* Fatal error: energy nonconservation:',
     &                '*********'
            WRITE(0,7899)DELTA_E
            WRITE(0,*)'**********************************************',
     &                '*********'
            
!            STOP
         ENDIF

      ENDIF

!*** diagnostic
!      write(0,*)'in vsg_td_emission, a=',a
!      write(0,*)'about to return to vsg_stat_therm'
!***

      RETURN
 6511 FORMAT(A3,' = composition',/,
     &       4I5,' = NC, NH, ND, ZG')
 6512 FORMAT(A3,' = composition',/,/)
 6521 FORMAT(A6,' = radiation field',/,/)
 6522 FORMAT(A6,' = radiation field',/,
     &       3E10.3,' = T_star (K), L/Lsol , distance (cm)')
 7601 FORMAT('try Umin/hc, Umax/hc=',1P2e12.5,' cm-1')
 7610 FORMAT('Pmax=',1PE12.5,' for J=',I3,' U/hc=',0PF11.2,' cm-1',' T=',
     &       0PF8.2,' K')
 7620 FORMAT('Umax/hc=',1PE12.5,' cm-1')
 7630 FORMAT('set UMAXHI=',1PE12.5,' cm-1')
 7650 FORMAT(' -- reduce UMAX -- , JCUT=',I4,
     &       ' new UMAX/hc=',F11.2,' cm-1')
 7670 FORMAT('increase umax: ',
     &       'umaxhi=',1PE12.5,' set umax=',0PF5.3,'*umax=',1PE12.5)
 7671 FORMAT('increase umax: ',
     &       'umaxhi=',1PE12.5,' set umax=1.1*umax=',1PE12.5)
 7680 FORMAT('increase umax: ',
     &       'umaxhi=',1PE12.5,' set umax=',1PE12.5)
 7690 FORMAT('increase umin: ',
     &       'JCUT=',I4,' uminlo=',1PE12.5,' set new umin=',1PE12.5)
 7694 FORMAT('decrease umin: ',
     &       'uminlo=',1PE12.5,' set new umin=0.95*uminlo=',1PE12.5)
 7695 FORMAT('decrease umin: ',
     &       'uminlo=',1PE12.5,' set new umin=0.9*uminlo=',1PE12.5)
 7696 FORMAT('decrease umin: ',
     &       'uminlo=',1PE12.5,' set new umin=0.8*uminlo=',1PE12.5)
 7700 FORMAT('Have iterated ',I5,' times, with error=',1PE10.3,/,
     &       'Matrix has been solved with TOL=',1PE10.3)
 7701 FORMAT('Bicg: sum_P=',1PE12.5)
 7710 FORMAT(A,/,A,/,A,/,A,/,A,/,A,/,A,/,A,/,A,/,
     &       '   lambda   nu*P_nu   C_abs',/,
     &       '   (um)    (erg s-1)  (cm2)')
 7750 FORMAT(F9.4,2(1PD12.4))
 7760 FORMAT('Prad=',1PE10.3,' erg s-1')
 7770 FORMAT('solved matrix: total emission Prad=',1PE10.3,' erg s-1')
 7771 FORMAT('single T=',F8.3,' K: Prad=',1PE10.3,' erg s-1')
 7780 FORMAT('size=',F8.2,
     &       ' A : powers computed from matrix transitions:',/,
     &       '  P_abs=',1PE10.3,' P_em=',1PE10.3,
     &       ' P_em(NSTATE)=',1PE10.3,' erg s-1')
 7800 FORMAT('ISRFabs=',1PE10.3,' frac error=',1PE10.2,/,50('='))
 7899 FORMAT(' ********* DELTA_E=',0PF12.5)
 7900 FORMAT('F=',I3,' P(F)=',1PE10.3,' norm.error=',1PE10.3)
 7901 FORMAT('energy check: E_em=',1PE10.3,' E_abs=',1PE10.3,
     &       ' E_em/E_abs-1=',1PE10.2)
 9876 format('PSTATE(',I3,')/PMAX=',1PE10.3,' U/hc=',0PF11.2,
     &       ' cm-1 T=',0PF8.2,' K')
      END
