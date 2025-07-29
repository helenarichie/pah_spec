      SUBROUTINE QCOMP(A,TGR,WAVE,ICOMP,QABS,QEXT,QSCA,G,G2,INDOUT)
      IMPLICIT NONE
c========================= qcomp v11 ===============================
C Parameters:

      INTEGER MXNANG,NCMAX
      PARAMETER(MXNANG=1000,NCMAX=3)

C Arguments:

      INTEGER ICOMP,INDOUT
      REAL A,G,G2,TGR,QABS,QEXT,QSCA,WAVE

C Common:
      DOUBLE PRECISION FIP
      COMMON/FIP_COM/FIP

C Local variables:

      INTEGER MOMDIM,MXANG,MXANGBH,NANGLMX

      PARAMETER(MOMDIM=3)
      PARAMETER(MXANG=MXNANG-1,MXANGBH=2*MXNANG-1)
      PARAMETER(NANGLMX=361)

      LOGICAL ANYANG,MAGNETIC,PERFCT,PRNT(2)
      INTEGER 
     &   ICOMPPA,ICOMPPE,ICOMPT,IDIEL_CPA,IDIEL_CPE,IPOLZN,J,JJ,
     &   N_C,N_H,NANG,NC,NMOM,NANGL,NTHETAS,ZPAH
      INTEGER
     &   IDIEL(1:NCMAX)
      REAL 
     &   AT,AT2,ENIM,ENRE,ENRE1,EPS11,EPS2,FPAH,G2PA,G2PE,GPA,GPE,M1ADT,
     &   PI,QABSPA,QABSPAH,QABSPE,QEXTPA,QEXTPE,QSCAPA,QSCAPE,
     &   X,XDIPL,XYMT,Y,Z
      REAL
     &   FC(1:NCMAX),
     &   GG(1:NCMAX),
     &   GG2(1:NCMAX),
     &   QQABS(1:NCMAX),
     &   QQEXT(1:NCMAX),
     &   QQSCA(1:NCMAX)
      DOUBLE PRECISION
     &   ALPHA,BETA,CABS_PAH,CHI0_DW,CHI0_PERP,CS,
     &   D2EDTHETA2,DPI,DTOH,DXX,DZZ,FAC,FGMIN,
     &   G2D,GD,GAMMA_A,GAMMA_B,GQSC,H_KAZ,H_KBZ,LAMBDA,LXX,LZZ,MIMCUT,
     &   M_AZ,M_BZ,M_S,N_AB,
     &   OMEGA,OMEGA_0,OMEGA_0A,OMEGA_0B,
     &   OMEGA_M,OMEGA_MA,OMEGA_MB,OMEGA0_DW,OMEGA0_PERP,
     &   QABSD,QEXTD,QSCAD,RADIUS,SN,SPIKE,SUM,TAU_DW,TAU_PERP,TGR_D,XX,
     &   WAVE_CM,WCABSMAGV,WGT
      DOUBLE PRECISION
     &   DQDO(0:MXANG),
     &   PMOM(0:MOMDIM,4),
     &   POL(1:NANGLMX),
     &   S11(1:NANGLMX),
     &   S12(1:NANGLMX),
     &   S33(1:NANGLMX),
     &   S34(1:NANGLMX),
     &   THETA(1:NANGLMX),
     &   THETAS(0:MXANG),
     &   XMU(1:MXANG)
      COMPLEX CXINDX,CXMU
      DOUBLE COMPLEX 
     &   DCXALPHAV_EDDY,DCXSBACK,DCXCHI_MINUS,DCXCHI_PLUS,DCXEPS,
     &   DCXSFORW,DCXI,DCXMU,DCXOMEGA_0AP,DCXOMEGA_0BP,DCXOMEGAPRIME,
     &   DCXPHI_EDDY,DCXPSI,DCXREF
      DOUBLE COMPLEX
     &   S1DCX(1:MXANG),
     &   S2DCX(1:MXANG),
     &   TBACK(2),
     &   TFORW(2)

c if decide to enable the section calling BHMIE, 
c need to declare following variables:

      REAL QBACK
      COMPLEX
     &   S1(MXANGBH),
     &   S2(MXANGBH)

c if decide to enable code calculating integrals over differential
c scattering cross section, need to declare following variables:
c      INTEGER JA
c      DOUBLE PRECISION DSUM,G1,S11,SUM
c      DOUBLE PRECISION
c     &   WGT(1:MXANG),

      EXTERNAL INDEX
C***********************************************************************
C
C Subroutine QCOMP computes optical properties of dust grains
C using:
C       subroutine INDEX to obtain refractive index

C       subroutine EMDIPL to compute optical properties of spheres small
C                  compared to the wavelength using electric and magnetic
C                  dipole cross sections (x < XDIPL)

C       subroutine MIEV0 to compute optical properties using Mie theory,
C                  for x + |m-1|x < XYMT

C       subroutine ADT to compute optical properties using anomalous
C                  diffraction theory when
C                  x >> 1 and |m-1| << M1ADT

C       subroutine GEOMOPT to compute optical properties of spheres using
C                  geometric optics, when none of above are possible

C Note: Rayleigh-Gans approximation is not employed.
C
C Input: A    = grain radius (micron)
C               for PAH/graphitic material (ICOMP=11-16) this is the
C               "graphite-equivalent" radius given by
C                  A = (N_C/470.)**(1./3.)
C               where N_C is the number of C atoms.  This relation is 
C               for an assumed density rho=2.24 g cm-3 for graphite, 
C               amu=1.6605e-24 g, and atomic wt = 12.011
C        TGR  = grain temperature (degK)
C        WAVE = wavelength (micron)
C        ICOMP=-3 for 2007 graphite, 1/3-2/3 approximation
C                     (2-component free-electron model for E||c with
C                      sigma_dc=62.5 mho/cm)
C              -2 for D03 graphite, 1/3-2/3 approximation
C                     (sigma_dc=62.5 mho/cm for E||c)
C              -1 for graphite, 1/3-2/3 approximation, eps from DL84
C               1 for graphite, E || c , epsilon from DL84
c                     (sigma_dc=29.1 mho/cm)
C               2 for graphite, E perp c, epsilon from DL84
C               3 for astronomical silicate, epsilon from DL84
C               4 for diamond, type IB, N/C=0.003
C               5 for meteoritic diamond, rho=2.3 g cm-3
C               6 for alpha-SiC
C               7 for astronomical silicate, smoothed UV
C               8 for new astronomical silicate
C               9 for new graphite, E || c axis, with
C                               sigma_dc=62.5 mho/cm
C              10 for crystalline H2O ice
C              11 for neutral-PAH/D03 graphitic material with D/H=0
C                               sigma_dc=62.5 mho/cm for E||c
C              12 for ionized-PAH/D03 graphitic material with D/H=0
C                               sigma_dc=62.5 mho/cm for E||c
C              13 for neutral-PAH/D03 graphitic material with D/H=0.25
C                               sigma_dc=62.5 mho/cm for E||c
C              14 for ionized-PAH/D03 graphitic material with D/H=0.25
C                               sigma_dc=62.5 mho/cm for E||c
C              15 for neutral-PAH/D03 graphitic material with D/H=0.5
C                               sigma_dc=62.5 mho/cm for E||c
C              16 for ionized-PAH/D03 graphitic material with D/H=0.5
C                               sigma_dc=62.5 mho/cm for E||c
C              17 for D03 astronomical silicate
C              18 for D03 graphite E || c, with sigma_dc=62.5 mho/cm
C              19 for D03 graphite E perpendicular to c axis
C              20 for D03 astronomical silicate with porosity(a)
C              21 for D03 graphite, E parallel to c, with porosity(a)
C              22 for D03 graphite, E perp. to c, with porosity(a)
C              23 for neutral PAH/D03 graphite (randomly-oriented), with
c                     porosity(a)
C              24 for ionized PAH/D03 graphite (randomly-oriented), with
C                     porosity(a)
C              25 for graphite, E para c, Djurisic & Li 99, opt data
C              26 for graphite, E para c, Djurisic & Li 99, EELS data
C              27 for graphite, E perp c, Djurisic & Li 99
C              28 for cellulose pyrolyzed at 800C (Jaeger etal 1998)
C                     for E.le.6.972eV, extrapolated to submm 
C              29 for cellulose pyrolyzed at 600C (Jaeger etal 1998)
C                     for E.le.6.214eV, extrapolated to submm
C              30 for amorph.car. AC1 from Rouleau & Martin 1991
C              31 for amorph.car. BE1 from Rouleau & Martin 1991
C              32 for amorph.car. FC21PS from Rouleau & Martin 1991
C              33 for amorph.car. HAPS from Rouleau & Martin 1991
c              34 for amorph.car. BE from Zubko etal 1996
c              35 for amorph.car. ACAR from Zubko etal 1996
c              36 for amorph.car. ACH2 from Zubko etal 1996
c              37 for amorph.car. amcACH2 (ZBCM96 modified and extended)
c              38 for neutral PAH/amcACH2
c              39 for ionized PAH/amcACH2
c              40 for neutral PAH/graphite/amcACH2
c                     D03 graphite with sigma_dc(E||c)=62.5 mho/cm
c              41 for ionized PAH/graphite/amcACH2
c                     D03 graphite with sigma_dc(E||c)=62.5 mho/cm
c              42 for D03 graphite, E||c, with new two-component
c                     free-electron model with sigma_dc(E||c)=62.5 mho/cm
c              43 for neutral PAH/graphite with new two-component
c                     free-electron model for E||c, sigma_dc=62.5 mho/cm
c                     and D/H = 0
c              44 for ionized PAH/graphite with new two-component
c                     free-electron model for E||c, sigma_dc=62.5 mho/cm
c                     and D/H = 0
c              45 for neutral PAH/graphite with new two-component
c                     free-electron model for E||c, sigma_dc=62.5 mho/cm
c                     and D/(H+D) = 0.25
c              46 for ionized PAH/graphite with new two-component
c                     free-electron model for E||c, sigma_dc=62.5 mho/cm
c                     and D/(H+D) = 0.25
c              47 for neutral PAH/graphite with new two-component
c                     free-electron model for E||c, sigma_dc=62.5 mho/cm
c                     and D/(H+D) = 0.5
c              48 for ionized PAH/graphite with new two-component
c                     free-electron model for E||c, sigma_dc=62.5 mho/cm
c                     and D/(H+D) = 0.5
c              49 for neutral PAH *** with no graphite ***
c              50 for ionized PAH *** with no graphite ***
c              51 for metallic Fe sphere (pure Fe, T=293K), no mag response
c              52 for metallic Fe sphere (pure Fe, T=293K), with magresponse
c              53 for magnetite Fe3O4 at T < 100K, no magnetic response
c              54 for magnetite Fe3O4 at T < 100K, with magnetic response
c              55 for maghemite Fe2O3 at T < 250K, no magnetic response
c              56 for maghemite Fe2O3 at T < 250K, with magnetic response
c              57 for metallic Fe 1:1:1.5 spheroid, with magnetic response
c              58 for metallic Fe 1:1:2 spheroid, with magnetic response
c              59 for metallic Fe 1:1:3 spheroid, with magnetic response
c              60 for metallic Fe 1:1:4 spheroid, with magnetic response
c              61 for metallic Fe 1:1:5 spheroid, with magnetic response
c              62 for metallic Fe 1:1:10 spheroid, with magnetic response
c              63 for metallic Fe 1:1:inf spheroid, with magnetic response
c        INDOUT=0 for no output, 1 to output refr. index
c
c Output: QABS =grain absorption efficiency
C         QEXT =grain extinction efficiency
C         QSCA =grain scattering efficiency
C          G   =<cos(theta)>
C          G2  =<cos^2(theta)>
C
C Control parameters (set in DATA statement below):
C         XIDPL : use dipole theory for x < XDIPL
C         XYMT  : use Mie theory for x > XDIPL and |m-1|x + x < XYMT
C         M1ADT : use ADT if x > XDIPL 
C                            and cannot use Mie theory
C                            and |m-1| < M1ADT
C         FGMIN : "graphitic" contribution to ICOMP=11 or 12 in limit
C                 a -> 0
C         AT   : "transition" radius from PAH-like to graphitic
C                 for ICOMP=11 or 12
C         AT2  : "trnasition" radius where amc composition begins

C Model for graphite fraction for PAH-graphite model:
c         V_pah = (1.-fgmin)*(4*pi*a^3/3)                for a < a_t
c               = (1.-fgmin)*(4*pi*a_t^3/3)              for a > a_t
c         V_gra = fgmin*(4*pi*a^3/3)                     for a < a_t
c               = fgmin*(4*pi*a_t^3/3)+(4*pi/3)*(a^3-a_t^3)  for a > a_t
c or
c         f_pah = (1.-fgmin)                             for a < a_t
c               = (1.-fgmin)*(a_t/a)^3                   for a > a_t
c         f_gra = fgmin                                  for a < a_t
c               = fgmin*(a_t/a)^3 + [1-(a_t/a)^3]        for a > a_t
c                    = 1 - (1-fgmin)*(a_t/a)^3

c Model for PAH-graphite-amc model
c         V_pah = (1.-fgmin)*(4*pi*a^3/3)                for a < a_t
c               = (1.-fgmin)*(4*pi*a_t^3/3)              for a > a_t
c         V_gra = fgmin*(4*pi*a^3/3)                     for a < a_t
c               = fgmin*(4*pi*a_t^3/3)+(4*pi/3)*(a^3-a_t^3) 
c                                                        for a_t < a < a_t2
c               = fgmin*(4*pi*a_t^3/3)+(4*pi/3)*(a_t2^3-a_t^3) 
c                                                        for a > a_t2
c         V_amc = 0                                      for a < a_t2
c               = (4*pi/3)*(a^3-a_t2^3)                  for a > a_t2
c or
c         f_pah = (1.-fgmin)                             for a < a_t
c               = (1.-fgmin)*(a_t/a)^3                   for a > a_t
c         f_gra = fgmin                                  for a < a_t
c               = fgmin*(a_t/a)^3+1-(a_t/a)^3            for a_t < a < a_t2
c                     = 1 - (1-fgmin)*(a_t/a)^3
c               = fgmin*(a_t/a)^3+[(a_t2)^3-(a_t)^3]/a^3 for a > a_t2
c                                                  
c         f_amc = 0                                      for a < a_t2
c               = 1-(a_t2/a)^3                           for a > a_t2
C***
C B.T.Draine, Princeton Univ. Obs., 89.1.25
C History:
C 91.07.16 (BTD) Modified to use BHMIE instead of MIETHE
C 91.07.16 (BTD) Modified to use GEOMOPT for large values of x
C 91.07.17 (BTD) Completed modifications.
C 91.07.30 (BTD) Change to |m|x as criterion for switching from Mie
C                theory to geometric optics:
C 91.07.30 (BTD) Add use of EMDIPL when x << 1.
C 91.07.31 (BTD) Further modifications
C 91.12.13 (BTD) Added handling of SiC
C 91.12.23 (BTD) Modified to use Rayleigh-Gans theory when
C                x > XDIPL but |m-1|x < ZRG
C 92.01.07 (BTD) Change criteria for choosing between
C                (1) dipole approximation
C                (2) Rayleigh-Gans theory
C                (3) Mie series expansion
C                (4) geometric optics
C                introduce variable NANG
C                set NANG=2 for calls to BHMIE
C 92.01.15 (BTD) changed XDIPL from .01 to .001
C                changed ZRG from .01 to .001
C 99.02.17 (BTD) added SAVE statement to make f77 compliant, and
C                g77 compatible
C 99.02.18 (BTD) modified to allow "astronomical silicate, smoothed UV"
C                as option
C 99.10.06 (BTD) modified to allow "PAH/graphitic material" as option
C                using subroutine PAH_crs_sct.f created by Aigen Li
C 99.10.08 (BTD) modified for ICOMP=11 for neutral PAH properties,
C                and ICOMP=12 for "ionized" PAH properties
C 99.10.25 (BTD) corrected error in calculation of N_H
C                set carbon density to 1.5 g cm-3 rather than 1.56
C 99.11.01 (BTD) groomed to please FTNCHEK
C 99.12.22 (BTD) modified to use fraction FGMIN of graphite for ICOMP=11
C                or 12 when radius A < AT
C 00.03.09 (BTD) modified to recognize new grain types ICOMP=13-16
C                for neutral,ionized PAHs with D/H=0.5 and 1
C                new variable introduced: DTOH = D/H ratio
C                DTOH is now argument of PAH_CRS_SCT
C 00.03.20 (BTD) reduced ZRG from 0.001 to 0.0005 to suppress small
C                discontinuity in computed cross section for very small
C                silicate grains in near IR (e.g., a=0.00035 and 
C                lambda=1.452um)
C 00.04.16 (BTD) changed definition of "radius" for PAH/graphitic
C                grains defined to be radius of grain with carbon
C                density of graphite
C 00.07.09 (BTD) Modified to use new version of pah_crs_sct provided
C                by Aigen Li for ICOMP=11-16:
C                for 1/lambda = 17.25 um-1 designed to reproduce
C                absorption/C for graphite.  For 1/lambda > 17.25 um-1
C                we use graphite for ICOMP=11-16
C 00.07.17 (BTD) Now can do ICOMP=8 (new astrosilicate)
C 00.09.13 (BTD) Add options ICOMP=-2 and ICOMP=9 to experiment with 
C                altering the FIR behavior of graphite with E parallel c
C                Modified ICOMP=11-16 so they now use the modified
C                graphite properties.
C                Note that only difference between ICOMP=1 and 9
C                is free electron contribution calculated by subroutine
C                FREE
C 00.10.27 (BTD) Modified ICOMP=11-16 so they again use DL84 graphite
C                properties
C                set FGMIN=0.01
C 02.10.02 (BTD) Modified ICOMP=11-16 so that QSCA is always taken to
C                be the value calculated for graphite (formerly we were
C                taking QSCA=FGMIN*QSCA(graphite) for very small PAHs
C                which is certainly an underestimate.
C                Note that this change has little practical impact
C                because for very small PAHs the scattering cross
C                sections are in any case very small.
C 02.12.22 (BTD) Modified to support ICOMP=17-23, where large grains
C                are allowed to be porous.
C                Porosity factor is set in subroutine INDEX
C 03.01.01 (BTD) Modified to support ICOMP=31-33 (experimental indices)
C 03.01.02 (BTD) Rewrote to have common module dealing with 1/3-2/3
C                approximation for graphite.
C                Modify to recognize new definitions of ICOMP=20-26
C 03.04.24 (BTD) Modified to use MIEV0 and ADT
C                No longer use RAYGAN
C                Add output to device 0 when BHMIE and MIEV0 differ
c                by more than 0.1% for Qext or Qsca
c 03.06.16 (BTD) Added calculation of <cos^2(theta)>, which requires
c                first calculating scattering phase function, and then
c                integrating over it.
c 03.06.16 (BTD) Added code (now commented out) which in future could
c                be used to calculate integrals over scattering phase
c                function.  At this time we do not need this because
c                MIEV0 already has provision for accurate calculation
c                of integrals over Legendre polynomials, which can
c                be used to obtain <cos> and <cos^2>
c 03.09.09 (BTD) Modify to clarify selection of method
c                Change criterion for use of Mie theory from
c                x < XYMT and |m|x < XYMT
c                to x + |m-1|x < XYMT
c                Remove variable ZRG (previously used for Rayleigh-Gans)
c                No longer any role for Rayleigh Gans approximation.
c 03.09.23 (BTD) Add M1ADT to SAVE statement.
c 03.09.25 (BTD) Corrected bugs concerning calculation of QABS from
c                QEXT and QSCA.
c                Added new double precision variable QABSD for internal
c                use.
c 03.09.28 (BTD) Corrected typo in computation of QABS for carbonaceous
c                grains.
c 04.06.02 (BTD) Added code to ensure that G2PA and G2PE are defined
c                even when in Rayleigh limit or when ADT or geometric
c                optics approximation is employed.  Note that in 
c                case of geometric optics G2PA and G2PE are not
c                accurate.
c 05.09.20 (BTD) changed dimensioning of S1 and S2 so that QCOMP and
c                BHMIE will be compatible when same value is assigned
c                to MXNANG in each.  Note that S1 and S2 are not used
c                here except as elements in argument list of BHMIE
c 06.06.02 (BTD) create version v2
c                corrected bug (typo) in treatment of ICOMP=13
c                changed options ICOMP=11-13 to now use D03 graphite
c                instead of DL84 graphite 
c 06.06.05 (BTD) new deuteration options:
c                changed ICOMP=13,14 from D/H=0.5 to D/H=0.25
c                changed ICOMP=15,16 from D/H=1 (100% deuterated)
c                                    to D/H=0.5 (50% deuterated)
c 06.08.24 (BTD) added new option:
c                ICOMP=27, randomly-oriented D03 graphite (no PAH)
c                now check that ICOMP has valid value
c 06.11.18 (BTD) add new options 
c                36-39 = AC1,BE1,FC21PS,HAPS from Rouleau & Martin 1991
c                40-42 = BE,ACAR,ACH2 from Zubko etal 1996
c 07.11.16 (BTD) add new option
c                43 = amcACH2 (ACH2 from ZBCM96 modified and extended)
c 07.11.26 (BTD) v3
c                * major rewrite to reorganize calculation
c                * add new options ICOMP=44,45,46,47
c 07.12.21 (BTD) v4
c                * delete old ICOMP=20,21
c                * renumber options
c                * use new numbering of IDIEL options consistent with
c                  new version (v4) of subroutine index
c 07.12.27 (BTD) * corrected major error -- had been computing PAH
c                  contribution when wave < 579.71A, whereas
c                  should compute this only when WAVE > 579.71A
c 08.01.04 (BTD) v5
c                * add new options:
c                  49 = pure PAH neutral (no graphite free-electron contrib.)
c                  50 = pure PAH ion (no graphite free-electron contrib.)
c                * restructured computation
c                * correct error in computation of PAHs in X-rays:
c                  had IDIEL_CPA=12,IDIEL_CPE=30
c                  should be IDIEL_CPA=30,IDIEL_CPE=13
c                  [this has no practical impact, since 2 dielectric functions
c                  are nearly identical at short wavelenghts]
c 08.07.21 (BTD) v6
c                * add COMMON/FIP_COM/FIP to communicate fraction of
c                  pah absorption cross section contributed by in-plane
c                  vibrational modes.  FIP is initialized to 2/3 here,
c                  but will then be reset by qpah_crs_sct_v15.f
c 11.02.23 (BTD) * modify to handle ICOMP=51, IDIEL=31
c 11.03.01 (BTD) * changed XDIPL from 1.e-3 to 1.e-4
c 11.03.05 (BTD) v7
c                * for Fe (ICOMP=51) add effect of magnetic permeability
c                  on eddy current dissipation
c                * for Fe (ICOMP=51) add magnetic dipole absorption
c                  due to Im(mu)
c                * fixed bugs
c                * changed xdipl from 1e-4 to 1e-2
c 11.03.12 (BTD) * incorporate call to subroutine MAGMIE to
c                  evaluate Qext, Qsca, Qabs for magnetic spheres
c                  [MAGMIE calls routines MAGSPH and MSPHSX written
c                   by M.E. Milham and recoded by Thomas Wriedt]
c                * added code to evaluate <cos^2> for scattering
c                  by magnetic sphere
c                * changed CXINDEX back to COMPLEX (from DOUBLE COMPLEX)
c 11.05.24 (BTD) * turn off magnetic absorption calculation for
c                  Fe by setting MAGNETIC=.FALSE.
c                * define DCXEPS and DCXMU as double complex
c 11.05.25 (BTD) * metallic Fe: modify to compute isotropic mu
c                  and use with magmie
c 11.06.15 (BTD) v8
c                * create new composition options:
c                  ICOMP=51 = Fe without magnetic response
c                        52      with
c                  ICOMP=53 = Fe3O4 without magnetic response
c                        54         with          
c 11.07.03 (BTD) * further edits to modify treatment of ferrimagnetic
c                  Fe3O4
c 11.07.11 (BTD) * changed N_AB from 100 to 10800 for Fe3O4
c 11.09.15 (BTD) v9
c                * create new composition options
c                  ICOMP=55 = maghemite Fe2O3 without magnetic response
c                        56 = maghemite Fe2O3 with magnetic response
c 11.10.12 (BTD) * set Gilbert damping parameter alpha_G = 0.2 for 
c                  Fe, Fe3O4, and Fe2O3
c 11.12.13 (BTD) * eliminate calculation of isotropized magnetic
c                  permeability
c                * for all materials to date, set MAGNETIC=.FALSE.
c                  so that Mie theory calculation is done with mu=1
c                * calculate magnetic contribution to absorption as
c                  a perturbation that can be added to Mie theory
c                  result
c 11.12.15 (BTD) * modify to assume 1:1:2 prolate spheroid when
c                  computing magnetization response for metallic Fe
c 12.05.07 (BTD) v10
c                * modify to assume sphere for ICOMP=52 magnetic Fe
c                * ICOMP=57: magnetic Fe with 1:1:1.5 prolate spheroid
c                * ICOMP=58: magnetic Fe with 1:1:2 prolate spheroid
c                * ICOMP=59: magnetic Fe with 1:1:3 prolate spheroid
c                * ICOMP=60: magnetic Fe with 1:1:4 prolate spheroid
c                * ICOMP=61: magnetic Fe with 1:1:5 prolate spheroid
c                * ICOMP=62: magnetic Fe with 1:1:10 prolate spheroid
c                * ICOMP=63: magnetic Fe with 1:1:inf prolate spheroid
c                * NB: for ICOMP=57-63 only allowance for nonsphericity
c                  is for calculation of ferrormagnetic resonance
c                  frequency omega_0, ranging from 5.9 GHz (for 1:1:1.5)
c                  to 32.4 GHz (for 1:1:inf)
c 12.12.08 (BTD) v11
c               * modify to include suppression of magnetic absorption
c                 due to eddy current shielding
c               * use routine ALPHA_EDDY to compute magnetic dipole
c                 polarizability due to eddy currents, and use this
c                 to estimate eddy current shielding correction
c                 phi_eddy
c               * add artificial limiter to keep suppress magnetic
c                 contribution when a/lambda > 0.001
c                 e.g., lambda=  1cm, a > 10um
c                                1mm, a > 1um
c                              100um, a > 0.1um
c                 under these circumstances
c                 (1) the magnetic calculation is not correct for
c                     multidomain grains
c                 (2) eddy current absorption completely dominates
c                     for Fe grains
! 18.07.06 (BTD) * /TABCOM/ -> /TABCOM2/
!                  TABLE -> TABLE2
!                  NTABMX=2000 -> NTABMX=10000
c end history
c***********************************************************************
      SAVE AT,AT2,DCXI,FGMIN,M1ADT,NANG,PI,XDIPL,XYMT

c Parameters for transition from PAH to graphite for ICOMP=11-16:
c A_T = 0.0050 micron  = "transition" radius
c FGMIN = 0.01         = minimum "graphitic" contribution

      DATA AT/0.0050D0/,AT2/0.0100D0/,FGMIN/0.01D0/
      DATA DCXI/(0.D0,1.D0)/,NANG/2/

C Parameters for selection of computational method:
C use dipole approx for x < XDIPL
C otherwise, use Mie Theory if x + |m-1|x < XYMT
C otherwise, use ADT if |m-1| < M1ADT

C M1ADT = max|m-1| for use of ADT
C XDIPL = max(x) for use of dipole approximation
C XYMT  = max(1+|m-1|x) for use of Mie theory

      DATA M1ADT/0.01/,XDIPL/0.01/,XYMT/20000./
      PI=4.*ATAN(1.E0)
      DPI=4.D0*ATAN(1.D0)
c=======================================================================
c*** diagnostic
c      write(0,*)'entered qcomp with a=',a
c      write(0,*)'                 tgr=',tgr
c      write(0,*)'                wave=',wave
c      write(0,*)'               icomp=',icomp
c      write(0,*)'              indout=',indout
c***
c Check that ICOMP is valid:

      IF(ICOMP.LT.-3.OR.ICOMP.GT.63)THEN
         WRITE(0,*)'qcomp_v11 ckpt 1:'
         WRITE(0,*)'   FATAL ERROR: QCOMP called with invalid ICOMP=',
     &             ICOMP
         STOP
      ENDIF

c initialize magnetic permeability to unity as default

      MAGNETIC=.FALSE.
      CXMU=(1.,0.)
      DCXMU=(1.D0,0.D0)

c initialize magnetic contribution to absorption to zero
c WCABSMAGV = lambda * C_abs,m / V

      WCABSMAGV=0.D0

C Compute scattering parameter X=2*pi*a/lambda

      X=2.*PI*A/WAVE

c=======================================================================
c set various input variables required by subroutine MIEV0:

      PERFCT=.FALSE.
      MIMCUT=0.D0
      ANYANG=.TRUE.

      NMOM=3

      IPOLZN=0
      PRNT(1)=.FALSE.
      PRNT(2)=.FALSE.

      NANG=1
      XMU(1)=0.

c***
c If following code is enabled,
c MIEV0 will calculate scattering in NANG directions, equally spaced
c in cos(theta) from 1 to -1
c This will allow direct calculation of integrals over differential
c scattering cross section.
c
c choose 
c      NANG=INT(MAX(20.,63.*A/WAVE))
c      NANG=MIN(NANG,MXNANG)
c
c have NANG be even to that we can use Simpson's rule
c
c      IF(NANG.LT.MXNANG)NANG=2*((NANG+1)/2)+1
c
c      DO JA=1,NANG
c         XMU(JA)=1.D0-DBLE(2*(JA-1))/DBLE(NANG-1)
c      ENDDO
c
c establish weights for Simpson's rule integration
c we assume that NANG > 2
c
c      WGT(1)=1.
c      WGT(2)=4.
c      DO JA=3,NANG-2,2
c         WGT(JA)=2.
c         WGT(JA+1)=4.
c      ENDDO
c      WGT(NANG)=1.
c***

c=======================================================================
c set quantities required by subroutine ADT:

      NTHETAS=0
      THETAS(0)=0.

c set certain defaults, which may later be changed
c Note: when pah_crs_sct_v15.f is called, it will reset the value of FIP
c       to the appropriate value for PAH neutral or ion.

      FPAH=0.
      DTOH=0.
      NC=1
      FC(1)=1.D0
      FIP=2.D0/3.D0

c=======================================================================
      IF(ICOMP.EQ.-3)THEN

c -3: graphite, 1/3-2/3 approximation
c     D03 graphite, with sigma_dc(E||c)=62.5 mho/cm
c                   with new 2 component free electron model for E||c
c                   to suppress 33um feature seen for IDIEL=12
         NC=2
         FC(1)=1./3.
         FC(2)=1.-FC(1)
         IDIEL(1)=30
         IDIEL(2)=13

c=======================================================================
      ELSEIF(ICOMP.EQ.-2)THEN

c -2: graphite, 1/3-2/3 approximation
c     D03 graphite, with sigma_dc(E||c)=62.5 mho/cm

         NC=2
         FC(1)=1./3.
         FC(2)=1.-FC(1)
         IDIEL(1)=12
         IDIEL(2)=13

c=======================================================================
      ELSEIF(ICOMP.EQ.-1)THEN

c -1: DL84 graphite, 1/3-2/3 approximation 
c                    sigma_dc(E||c)=29.1 mho/cm

         NC=2
         FC(1)=1./3.
         FC(2)=1.-FC(1)
         IDIEL(1)=1
         IDIEL(2)=2

c=======================================================================
      ELSEIF(ICOMP.EQ.3)THEN
c 3: DL84 astronomical silicate
         IDIEL(1)=3
c=======================================================================
      ELSEIF(ICOMP.EQ.4)THEN
c 4: diamond, type IB, N/C=0.003
         IDIEL(1)=4
c=======================================================================
      ELSEIF(ICOMP.EQ.5)THEN
c 5: meteoritic diamond, rho=2.3 g cm-3
         IDIEL(1)=5
c=======================================================================
      ELSEIF(ICOMP.EQ.6)THEN
c 6: alpha-SiC
         IDIEL(1)=6
c=======================================================================
      ELSEIF(ICOMP.EQ.7)THEN
c 7: astronomical silicate, smoothed UV
         IDIEL(1)=7
c=======================================================================
      ELSEIF(ICOMP.EQ.8)THEN
c 8: new astronomical silicate
         IDIEL(1)=8
c=======================================================================
      ELSEIF(ICOMP.EQ.9)THEN
c 9: DL84 graphite, E parallel to c axis, but
c    with sigma_dc(E||c)=62.5 mho/cm instead of 29.1 mho/cm
         IDIEL(1)=9
c=======================================================================
      ELSEIF(ICOMP.EQ.10)THEN
c 10: crystalline H2O ice
         IDIEL(1)=10
c=======================================================================
      ELSEIF(ICOMP.EQ.11)THEN
c 11: neutral-PAH/D03 graphitic material with D/H=0
c     with sigma_dc(E||c)=62.5 mho/cm
         ZPAH=0
         DTOH=0.
c for wavelengths < 579.71A, assume pure graphite:
         IF(WAVE.GT.0.057971)THEN
            IF(A.LE.AT)THEN
               FPAH=REAL(1.-FGMIN)
            ELSE
               FPAH=REAL(1.-FGMIN)*(AT/A)**3
            ENDIF
         ENDIF
         NC=2
         FC(1)=(1.-FPAH)/3.
         FC(2)=2.*FC(1)
         IDIEL(1)=12
         IDIEL(2)=13
c=======================================================================
      ELSEIF(ICOMP.EQ.12)THEN
c 12: for ionized-PAH/D03 graphitic material with D/H=0
c     with sigma_dc(E||c)=62.5 mho/cm
         ZPAH=1
         DTOH=0.
c for wavelengths < 579.71A, assume pure graphite:
         IF(WAVE.GT.0.057971)THEN
            IF(A.LE.AT)THEN
               FPAH=REAL(1.-FGMIN)
            ELSE
               FPAH=REAL(1.-FGMIN)*(AT/A)**3
            ENDIF
         ENDIF
         NC=2
         FC(1)=(1.-FPAH)/3.
         FC(2)=2.*FC(1)
         IDIEL(1)=12
         IDIEL(2)=13
c=======================================================================
      ELSEIF(ICOMP.EQ.13)THEN
c 13: neutral-PAH/D03 graphitic material with D/H=0.25
         ZPAH=0
         DTOH=0.25
c for wavelengths < 579.71A, assume pure graphite:
         IF(WAVE.GT.0.057971)THEN
            IF(A.LE.AT)THEN
               FPAH=REAL(1.-FGMIN)
            ELSE
               FPAH=REAL(1.-FGMIN)*(AT/A)**3
            ENDIF
         ENDIF
         NC=2
         FC(1)=(1.-FPAH)/3.
         FC(2)=2.*FC(1)
         IDIEL(1)=12
         IDIEL(2)=13
c=======================================================================
      ELSEIF(ICOMP.EQ.14)THEN
c 14: ionized-PAH/D03 graphitic material with D/H=0.25
         ZPAH=1
         DTOH=0.25
c for wavelengths < 579.71A, assume pure graphite:
         IF(WAVE.GT.0.057971)THEN
            IF(A.LE.AT)THEN
               FPAH=REAL(1.-FGMIN)
            ELSE
               FPAH=REAL(1.-FGMIN)*(AT/A)**3
            ENDIF
         ENDIF
         NC=2
         FC(1)=(1.-FPAH)/3.
         FC(2)=2.*FC(1)
         IDIEL(1)=12
         IDIEL(2)=13
c=======================================================================
      ELSEIF(ICOMP.EQ.15)THEN
c 15: neutral-PAH/D03 graphitic material with D/H=0.5
         ZPAH=0
         DTOH=0.50
c for wavelengths < 579.71A, assume pure graphite:
         IF(WAVE.GT.0.057971)THEN
            IF(A.LE.AT)THEN
               FPAH=REAL(1.-FGMIN)
            ELSE
               FPAH=REAL(1.-FGMIN)*(AT/A)**3
            ENDIF
         ENDIF
         NC=2
         FC(1)=(1.-FPAH)/3.
         FC(2)=2.*FC(1)
         IDIEL(1)=12
         IDIEL(2)=13
c=======================================================================
      ELSEIF(ICOMP.EQ.16)THEN
c 16: ionized-PAH/D03 graphitic material with D/H=0.5
         ZPAH=1
         DTOH=0.50
c for wavelengths < 579.71A, assume pure graphite:
        IF(WAVE.GT.0.057971)THEN
            IF(A.LE.AT)THEN
               FPAH=REAL(1.-FGMIN)
            ELSE
               FPAH=REAL(1.-FGMIN)*(AT/A)**3
            ENDIF
         ENDIF
         NC=2
         FC(1)=(1.-FPAH)/3.
         FC(2)=2.*FC(1)
         IDIEL(1)=12
         IDIEL(2)=13
c=======================================================================
      ELSEIF(ICOMP.EQ.17)THEN
c 17: D03 astronomical silicate
         IDIEL(1)=11
c=======================================================================
      ELSEIF(ICOMP.EQ.18)THEN
c 18: D03 graphite E parallel to c axis, with sigma_dc=62.5 mho/cm
         IDIEL(1)=12
c=======================================================================
      ELSEIF(ICOMP.EQ.19)THEN
c 19: D03 graphite E perpendicular to c axis
         IDIEL(1)=13
c=======================================================================
      ELSEIF(ICOMP.EQ.20)THEN
c 20: D03 astronomical silicate with porosity(a)
         IDIEL(1)=14
c=======================================================================
      ELSEIF(ICOMP.EQ.21)THEN
c 21: D03 graphite, E parallel to c, with porosity(a)
         IDIEL(1)=15
c=======================================================================
      ELSEIF(ICOMP.EQ.22)THEN
c 22: D03 graphite, E perp. to c, with porosity(a)
         IDIEL(1)=16
c=======================================================================
      ELSEIF(ICOMP.EQ.23)THEN
c 23: neutral PAH/D03 graphite (randomly-oriented), with porosity(a)
         ZPAH=0
c for wavelengths < 579.71A, assume pure graphite:
         IF(WAVE.GT.0.057971)THEN
            IF(A.LE.AT)THEN
               FPAH=REAL(1.-FGMIN)
            ELSE
               FPAH=REAL(1.-FGMIN)*(AT/A)**3
            ENDIF
         ENDIF
         NC=2
         FC(1)=(1.-FPAH)/3.
         FC(2)=2.*FC(1)
         IDIEL(1)=12
         IDIEL(2)=13
c=======================================================================
      ELSEIF(ICOMP.EQ.24)THEN
c 24: ionized PAH/D03 graphite (randomly-oriented), with porosity(a)
c                 sigma_dc(E||c)=62.5 mho/cm
         NC=2
         ZPAH=1
         DTOH=0.
c for wavelengths < 579.71A, assume pure graphite:
         IF(WAVE.GT.0.057971)THEN
            IF(A.LE.AT)THEN
               FPAH=REAL(1.-FGMIN)
            ELSE
               FPAH=REAL(1.-FGMIN)*(AT/A)**3
            ENDIF
         ENDIF
         FC(1)=(1.-FPAH)/3.
         FC(2)=2.*FC(1)
         IDIEL(1)=12
         IDIEL(2)=13
c=======================================================================
      ELSEIF(ICOMP.EQ.25)THEN
c 25: graphite, E para c, Djurisic & Li 99, opt data
         IDIEL(1)=17
c=======================================================================
      ELSEIF(ICOMP.EQ.26)THEN
c 26: graphite, E para c, Djurisic & Li 99, EELS data
         IDIEL(1)=18
c=======================================================================
      ELSEIF(ICOMP.EQ.27)THEN
c 27: graphite, E perp c, Djurisic & Li 99
         IDIEL(1)=19
c=======================================================================
      ELSEIF(ICOMP.EQ.28)THEN
c 28: cellulose pyrolyzed at 800C (Jaeger etal 1998) 
c     for E.le.6.972eV, extrapolated to submm 
         IDIEL(1)=20
c=======================================================================
      ELSEIF(ICOMP.EQ.29)THEN
c 29: cellulose pyrolyzed at 600C (Jaeger etal 1998)
c     for E.le.6.214eV, extrapolated to submm
         IDIEL(1)=21
c=======================================================================
      ELSEIF(ICOMP.EQ.30)THEN
c 30: amorph.car. AC1 from Rouleau & Martin 1991
         IDIEL(1)=22
c=======================================================================
      ELSEIF(ICOMP.EQ.31)THEN
c 31: amorph.car. BE1 from Rouleau & Martin 1991
         IDIEL(1)=23
c=======================================================================
      ELSEIF(ICOMP.EQ.32)THEN
c 32: amorph.car. FC21PS from Rouleau & Martin 1991
         IDIEL(1)=24
c=======================================================================
      ELSEIF(ICOMP.EQ.33)THEN
c 33: amorph.car. HAPS from Rouleau & Martin 1991
         IDIEL(1)=25
c=======================================================================
      ELSEIF(ICOMP.EQ.34)THEN
c 34: amorph.car. BE from Zubko etal 1996
         IDIEL(1)=26
c=======================================================================
      ELSEIF(ICOMP.EQ.35)THEN
c 35: amorph.car. ACAR from Zubko etal 1996 
         IDIEL(1)=27
c=======================================================================
      ELSEIF(ICOMP.EQ.36)THEN
c 36: amorph.car. ACH2 from Zubko etal 1996
         IDIEL(1)=28
c=======================================================================
      ELSEIF(ICOMP.EQ.37)THEN
c 37: amorph.car. amcACH2 (ZBCM96 modified and extended)
         IDIEL(1)=29
c=======================================================================
      ELSEIF(ICOMP.EQ.38)THEN
c 38: neutral PAH/amcACH2
         ZPAH=0
         IF(A.LE.AT2)THEN
            FPAH=REAL(1.-FGMIN)
         ELSE
            FPAH=REAL(1.-FGMIN)*(AT2/A)**3
         ENDIF
         FC(1)=1.-FPAH
         IDIEL(1)=29
      ELSEIF(ICOMP.EQ.39)THEN
c 39: ionized PAH/amcACH2
         ZPAH=1
         DTOH=0.
c for wavelengths < 579.71A, assume pure amcACH2
         IF(WAVE.GT.0.057971)THEN
            IF(A.LE.AT2)THEN
               FPAH=REAL(1.-FGMIN)
            ELSE
               FPAH=REAL(1.-FGMIN)*(AT2/A)**3
            ENDIF
         ENDIF
         FC(1)=1.-FPAH
         IDIEL(1)=29
c=======================================================================
      ELSEIF(ICOMP.EQ.40)THEN

c 40: neutral PAH/graphite/amcACH2
c     D03 graphite with two-component free-electron model
c                  with sigma_dc(E||c)=62.5 mho/cm
         NC=3
         ZPAH=1
         DTOH=0.
c for wavelengths < 579.71A, assume pure amcACH2
         IF(WAVE.GT.0.057971)THEN
            IF(A.LE.AT)THEN
               FPAH=REAL(1.-FGMIN)
            ELSE
               FPAH=REAL(1.-FGMIN)*(AT/A)**3
            ENDIF
         ENDIF
         IF(A.LE.AT2)THEN
            NC=2
            FC(3)=0.
         ELSE
            NC=2
            FC(3)=1-(AT2/A)**3
         ENDIF
         FC(1)=(1.-FPAH-FC(3))/3.
         FC(2)=2.*FC(1)
         IDIEL(1)=30
         IDIEL(2)=13
         IDIEL(3)=29
c=======================================================================
      ELSEIF(ICOMP.EQ.41)THEN
c 41: ionized PAH/graphite/amcACH2
         ZPAH=1
         DTOH=0.
c for wavelengths < 579.71A, assume pure amcACH2
         IF(WAVE.GT.0.057971)THEN
            IF(A.LE.AT)THEN
               FPAH=REAL(1.-FGMIN)
            ELSE
               FPAH=REAL(1.-FGMIN)*(AT/A)**3
            ENDIF
         ENDIF
         IF(A.LE.AT2)THEN
            NC=2
            FC(3)=0.
         ELSE
            NC=3
            FC(3)=1-(AT2/A)**3
         ENDIF
         FC(1)=(1.-FPAH-FC(3))/3.
         FC(2)=2.*FC(1)
         IDIEL(1)=30
         IDIEL(2)=13
         IDIEL(3)=29
c*** diagnostic
c         write(0,*)'set idiel(1-3)=18, 19, 29'
c         write(0,*)'fpah=',fpah
c         write(0,*)'fc(1)=',fc(1)
c         write(0,*)'fc(2)=',fc(2)
c         write(0,*)'fc(3)=',fc(3)
c***
c=======================================================================
      ELSEIF(ICOMP.EQ.42)THEN
c 42: D03 graphite, E||c, with new two-component
c         free-electron model with sigma_dc(E||c)=62.5 mho/cm
         IDIEL(1)=30
c=======================================================================
      ELSEIF(ICOMP.EQ.43)THEN
c 43: neutral PAH/D03 graphite (randomly-oriented) with new two-component
c         free-electron model with sigma_dc(E||c)=62.5 mho/cm
c         and D/H=0 
         NC=2
         ZPAH=0
         DTOH=0.
c for wavelengths < 579.71A, assume pure graphite
         IF(WAVE.GT.0.057971)THEN
            IF(A.LE.AT)THEN
               FPAH=REAL(1.-FGMIN)
            ELSE
               FPAH=REAL(1.-FGMIN)*(AT/A)**3
            ENDIF
         ENDIF
         FC(1)=(1.-FPAH)/3.
         FC(2)=2.*FC(1)
         IDIEL(1)=30
         IDIEL(2)=13
c=======================================================================
      ELSEIF(ICOMP.EQ.44)THEN
c 44: ionized PAH/D03 graphite (randomly-oriented) with new two-component
c         free-electron model with sigma_dc(E||c)=62.5 mho/cm
c         and D/H=0
         NC=2
         ZPAH=1
         DTOH=0.
c for wavelengths < 579.71A, assume pure graphite
         IF(WAVE.GT.0.057971)THEN
            IF(A.LE.AT)THEN
               FPAH=REAL(1.-FGMIN)
            ELSE
               FPAH=REAL(1.-FGMIN)*(AT/A)**3
            ENDIF
         ENDIF
         FC(1)=(1.-FPAH)/3.
         FC(2)=2.*FC(1)
         IDIEL(1)=30
         IDIEL(2)=13
c=======================================================================
      ELSEIF(ICOMP.EQ.45)THEN
c 45: neutral PAH/D03 graphite (randomly-oriented) with new two-component
c         free-electron model with sigma_dc(E||c)=62.5 mho/cm
c         and D/(H+D)=0.25
         NC=2
         ZPAH=0
         DTOH=0.25
c for wavelengths < 579.71A, assume pure graphite
         IF(WAVE.GT.0.057971)THEN
            IF(A.LE.AT)THEN
               FPAH=REAL(1.-FGMIN)
            ELSE
               FPAH=REAL(1.-FGMIN)*(AT/A)**3
            ENDIF
         ENDIF
         FC(1)=(1.-FPAH)/3.
         FC(2)=2.*FC(1)
         IDIEL(1)=30
         IDIEL(2)=13
c=======================================================================
      ELSEIF(ICOMP.EQ.46)THEN
c 46: ionized PAH/D03 graphite (randomly-oriented) with new two-component
c         free-electron model with sigma_dc(E||c)=62.5 mho/cm
c         and D/(H+D)=0.25
         NC=2
         ZPAH=1
         DTOH=0.25
c for wavelengths < 579.71A, assume pure graphite
         IF(WAVE.GT.0.057971)THEN
            IF(A.LE.AT)THEN
               FPAH=REAL(1.-FGMIN)
            ELSE
               FPAH=REAL(1.-FGMIN)*(AT/A)**3
            ENDIF
         ENDIF
         FC(1)=(1.-FPAH)/3.
         FC(2)=2.*FC(1)
         IDIEL(1)=30
         IDIEL(2)=13
c=======================================================================
      ELSEIF(ICOMP.EQ.47)THEN

c 47: neutral PAH/D03 graphite (randomly-oriented) with new two-component
c         free-electron model with sigma_dc(E||c)=62.5 mho/cm

         NC=2
         ZPAH=0
         DTOH=0.5

c for wavelengths < 579.71A, assume pure graphite

         IF(WAVE.GT.0.057971)THEN
            IF(A.LE.AT)THEN
               FPAH=REAL(1.-FGMIN)
            ELSE
               FPAH=REAL(1.-FGMIN)*(AT/A)**3
            ENDIF
         ENDIF
         FC(1)=(1.-FPAH)/3.
         FC(2)=2.*FC(1)
         IDIEL(1)=30
         IDIEL(2)=13
c=======================================================================
      ELSEIF(ICOMP.EQ.48)THEN

c 48: ionized PAH/D03 graphite (randomly-oriented) with new two-component
c         free-electron model with sigma_dc(E||c)=62.5 mho/cm

         NC=2
         ZPAH=1
         DTOH=0.5

c for wavelengths < 579.71A, assume pure graphite

         IF(WAVE.GT.0.057971)THEN
            IF(A.LE.AT)THEN
               FPAH=REAL(1.-FGMIN)
            ELSE
               FPAH=REAL(1.-FGMIN)*(AT/A)**3
            ENDIF
         ENDIF
         FC(1)=(1.-FPAH)/3.
         FC(2)=2.*FC(1)
         IDIEL(1)=30
         IDIEL(2)=13
c=======================================================================
      ELSEIF(ICOMP.EQ.49)THEN

c 49: neutral PAH (pure -- no graphite free-electron contribution)

         NC=0
         FPAH=1.
         ZPAH=0
         DTOH=0.

c=======================================================================
      ELSEIF(ICOMP.EQ.50)THEN

c 50: ionized PAH (pure -- no graphite free-electron contribution)
c     at high energies, use graphite, E perp c

         NC=0
         FPAH=1.
         ZPAH=1
         DTOH=0.

c=======================================================================
      ELSEIF(ICOMP.EQ.51)THEN

c 51: metallic Fe, without magnetic response

         IDIEL(1)=31

c=======================================================================
      ELSEIF(ICOMP.EQ.52.OR.(ICOMP.GE.57.AND.ICOMP.LE.63))THEN

c 52,57,58,59,60,61,62,63: metallic Fe, with magnetic response
c                          now include effects of eddy current shielding

         IDIEL(1)=31

c estimate magnetic permeability of metallic Fe

         OMEGA=2.*DPI*3.D14/WAVE

c domain wall contribution (see Draine & Lazarian 1999, ApJ 512, 740)

c         CHI0_DW=10.
c         TAU_DW=1.E-9
c         OMEGA0_DW=1.E10
c         CXDENOM=1.-(OMEGA/OMEGA0_DW)**2-DCXI*OMEGA*TAU_DW
c         DCXMU=1.+4.*PI*CHI0_DW/CXDENOM

c         CHI0_PERP=3.3
c         TAU_PERP=1.6E-11
c         OMEGA0_PERP=2.*PI*2.E10
c         CXDENOM=(1.-(OMEGA/OMEGA0_PERP)**2)-DCXI*OMEGA*TAU_PERP
c         DCXMU=DCXMU+4.D0*PI*CHI0_DW/CXDENOM

c for bcc Fe:
c    H_K = 548 G
c    4*pi*M_0 = 22020 Oe
c assume Gilbert damping parameter alpha = 0.2

         IF(ICOMP.EQ.52)THEN
            LZZ=1./3.
         ELSEIF(ICOMP.EQ.57)THEN
            LZZ=0.23298
         ELSEIF(ICOMP.EQ.58)THEN
            LZZ=0.17356
         ELSEIF(ICOMP.EQ.59)THEN
            LZZ=0.10871
         ELSEIF(ICOMP.EQ.60)THEN
            LZZ=0.07541
         ELSEIF(ICOMP.EQ.61)THEN
            LZZ=0.05582
         ELSEIF(ICOMP.EQ.62)THEN
            LZZ=0.02029
         ELSEIF(ICOMP.EQ.63)THEN
            LZZ=0.
         ENDIF
         DZZ=4.*PI*LZZ
         DXX=2.*PI*(1.-LZZ)
         ALPHA=0.2D0
         D2EDTHETA2=9.6D5
         M_S=22020./(4.*PI)

         H_KAZ=D2EDTHETA2/M_S
         GAMMA_A=-1.759D7
         OMEGA_0=-GAMMA_A*(H_KAZ-(DZZ-DXX)*M_S)
         OMEGA_M=-GAMMA_A*M_S
         DCXOMEGAPRIME=OMEGA_0-DCXI*ALPHA*OMEGA

c chi_\pm = omega_M/(omega_0' \pm omega)
c     psi = (1/2)(chi_+ + chi_-)
c         = (1/2)(2*omega_M * omega_0')/((omega_0')^2 - omega^2)
c         = omega_0' * omega_M / ((omega_0')^2 - omega^2)

         DCXPSI=DCXOMEGAPRIME*OMEGA_M/(DCXOMEGAPRIME**2-OMEGA**2)

c compute eddy current effects

         RADIUS=DBLE(A)
         TGR_D=DBLE(TGR)
         CALL ALPHA_EDDY(IDIEL(1),OMEGA,RADIUS,TGR_D,DCXALPHAV_EDDY)
         DCXPHI_EDDY=-(8.*DPI/3.)*DCXALPHAV_EDDY

c WCABSMAGV = lambda*C_abs,m/V
c where C_abs,m = magnetic contribution to absorption
c                 including shielding effects of eddy currents,
c                 but not including eddy current dissipation
c                 (eddy current dissipation will be calculated
c                  separately and added later)

         WCABSMAGV=((16./3.)*PI**2)*IMAG((1.D0-DCXPHI_EDDY)*DCXPSI)

c*** diagnostic
c         write(0,*)'qcomp_v11 ckpt 10'
c         write(0,*)'      omega_0=',omega_0
c         write(0,*)'      omega_M=',omega_M
c         write(0,*)'dcxomegaprime=',dcxomegaprime
c         write(0,*)'       dcxpsi=',dcxpsi
c         write(0,9100)WCABSMAGV
c***
c=======================================================================
      ELSEIF(ICOMP.EQ.53)THEN

c ICOMP=53: magnetite Fe3O4, no magnetic response

         IDIEL(1)=32

c=======================================================================
      ELSEIF(ICOMP.EQ.54)THEN

c ICOMP=54: magnetite Fe3O4, with magnetic response
c           assume Gilbert damping parameter alpha = 0.2

         IDIEL(1)=32
         MAGNETIC=.TRUE.

         OMEGA=2.*DPI*3.D14/WAVE

         GAMMA_A=-1.759D7
         ALPHA=0.2D0
         BETA=5.D0/9.D0
         D2EDTHETA2=29.2D5
         M_S=6400.D0/(4.D0*PI)
         N_AB=10800.D0
         GAMMA_B=GAMMA_A
         H_KAZ=((1.-BETA)/(1.+BETA))*D2EDTHETA2/M_S
         H_KBZ=-H_KAZ
         M_AZ=M_S/(1.-BETA)
         M_BZ=-BETA*M_AZ
         OMEGA_0A=-GAMMA_A*(H_KAZ-N_AB*M_BZ)
         OMEGA_0B=-GAMMA_B*(H_KBZ-N_AB*M_AZ)
         OMEGA_MA=-GAMMA_A*M_AZ
         OMEGA_MB=-GAMMA_B*M_BZ

         DCXOMEGA_0AP=OMEGA_0A-DCXI*ALPHA*OMEGA
         DCXOMEGA_0BP=OMEGA_0B+DCXI*ALPHA*OMEGA
         DCXCHI_PLUS=(OMEGA_MA*(DCXOMEGA_0BP-OMEGA)+
     &           OMEGA_MB*(DCXOMEGA_0AP-OMEGA)-
     &           2.*N_AB*OMEGA_MA*OMEGA_MB)/
     &           ((DCXOMEGA_0AP-OMEGA)*(DCXOMEGA_0BP-OMEGA)-
     &           N_AB**2*OMEGA_MA*OMEGA_MB)
         DCXCHI_MINUS=(OMEGA_MA*(DCXOMEGA_0BP+OMEGA)+
     &           OMEGA_MB*(DCXOMEGA_0AP+OMEGA)-
     &           2.*N_AB*OMEGA_MA*OMEGA_MB)/
     &           ((DCXOMEGA_0AP+OMEGA)*(DCXOMEGA_0BP+OMEGA)-
     &           N_AB**2*OMEGA_MA*OMEGA_MB)
         DCXPSI=(DCXCHI_PLUS+DCXCHI_MINUS)/2.D0

c the electrical conductivity of magnetite is low enough that
c we do not need to worry about eddy current shielding

         WCABSMAGV=((16./3.)*PI**2)*IMAG(DCXPSI)

c=======================================================================
      ELSEIF(ICOMP.EQ.55)THEN

c ICOMP=55: maghemite gamma-Fe2O3, no magnetic response

         IDIEL(1)=33

c=======================================================================
      ELSEIF(ICOMP.EQ.56)THEN

c ICOMP=56: maghemite gamma-Fe2O3, with magnetic response
c           assume Gilbert damping parameter alpha = 0.2

         IDIEL(1)=33
         MAGNETIC=.TRUE.
         
         OMEGA=2.*DPI*3.D14/WAVE

         GAMMA_A=-1.759D7
         ALPHA=0.2D0
         BETA=3.D0/5.D0
         D2EDTHETA2=10.9D5
         M_S=4890.D0/(4.*PI)
         N_AB=9280.D0

         GAMMA_B=GAMMA_A
         H_KAZ=((1.-BETA)/(1.+BETA))*D2EDTHETA2/M_S
         H_KBZ=-H_KAZ
         M_AZ=M_S/(1.-BETA)
         M_BZ=-BETA*M_AZ
         OMEGA_0A=-GAMMA_A*(H_KAZ-N_AB*M_BZ)
         OMEGA_0B=-GAMMA_B*(H_KBZ-N_AB*M_AZ)
         OMEGA_MA=-GAMMA_A*M_AZ
         OMEGA_MB=-GAMMA_B*M_BZ

         DCXOMEGA_0AP=OMEGA_0A-DCXI*ALPHA*OMEGA
         DCXOMEGA_0BP=OMEGA_0B+DCXI*ALPHA*OMEGA
         DCXCHI_PLUS=(OMEGA_MA*(DCXOMEGA_0BP-OMEGA)+
     &           OMEGA_MB*(DCXOMEGA_0AP-OMEGA)-
     &           2.*N_AB*OMEGA_MA*OMEGA_MB)/
     &           ((DCXOMEGA_0AP-OMEGA)*(DCXOMEGA_0BP-OMEGA)-
     &           N_AB**2*OMEGA_MA*OMEGA_MB)
         DCXCHI_MINUS=(OMEGA_MA*(DCXOMEGA_0BP+OMEGA)+
     &           OMEGA_MB*(DCXOMEGA_0AP+OMEGA)-
     &           2.*N_AB*OMEGA_MA*OMEGA_MB)/
     &           ((DCXOMEGA_0AP+OMEGA)*(DCXOMEGA_0BP+OMEGA)-
     &           N_AB**2*OMEGA_MA*OMEGA_MB)
         DCXPSI=(DCXCHI_PLUS+DCXCHI_MINUS)/2.D0

c electrical conductivity of maghemite is low, and we do not need to
c include effects of magnetic shielding by eddy currents

         WCABSMAGV=((16./3.)*PI**2)*IMAG(DCXPSI)
c=======================================================================
      ELSE
         WRITE(0,*)'qcomp_v11 ckpt 11: Fatal error: ICOMP=',ICOMP
         STOP
      ENDIF

c=======================================================================
c selection of method:

c if x < x_dipl  ---> EMDIPL

c else

c       if x < xymt and y=mx < xymt  ---> Mie theory

c       else if |m-1| < 0.01 ---> ADT

c       else ---> GEOMOPT

c       endif

c endif
c=======================================================================
c*** diagnostic
c      write(0,*)'qcomp_v11 ckpt 20: nc=',nc
c      write(0,*)'               fpah=',fpah
c      write(0,*)'               fc(1)=',fc(1)
c      write(0,*)'               fc(2)=',fc(2)
c      write(0,*)'               fc(3)=',fc(3)
c***
      QABS=0.
      QSCA=0.
      QEXT=0.
      G=0.
      G2=0.

c do we need to include PAH contribution?

      IF(FPAH.GT.0.)THEN

c Calculate PAH contribution to cross section

         IF(WAVE.GT.0.057971)THEN

C Use empirical model for PAH absorption
C Because N_C should be an integer, cannot consider PAH particles
C with more than 2**31 C atoms

            IF(A.LT.0.15)THEN
               N_C=NINT(4.70D11*DBLE(A**3))
            ELSE
               N_C=NINT(4.70D11*DBLE(0.15D0**3))
            ENDIF

C Note that coronene = C_24H_12 is pericondensed with H/C=0.5
C                      C_52H_18 is pericondensed with H/C=0.346
C                      C_96H_24 is pericondensed with H/C=0.25
C
C These are most extreme examples of pericondensation.
C We will assume H/C = 0.5 for N_C < 25
C                H/C = 0.5/(N_C/25)^{1/2} for 25 < N_C < 100
C                H/C = 0.25 for 100 < N_C 

            IF(N_C.LE.25.D0)THEN
               N_H=NINT(0.5D0*N_C)
            ELSEIF(N_C.GT.100.D0)THEN
               N_H=NINT(0.25D0*N_C)
            ELSE
               N_H=NINT(2.5D0*SQRT(DBLE(N_C)))
            ENDIF
            WAVE_CM=1.D-4*WAVE
            CALL PAH_CRS_SCT(N_C,N_H,DTOH,ZPAH,WAVE_CM,CABS_PAH)
            IF(A.LT.0.15)THEN
               QABS=FPAH*REAL(N_C*CABS_PAH/DBLE(1.E-8*PI*A**2))
            ELSE
               QABS=FPAH*4.70E19*A*REAL(CABS_PAH)/PI
            ENDIF
            QEXT=QABS

c*** diagnostic
c               write(0,*)'pah contribution to qabs=',qabs
c***
         ELSE

c at high energies, we assume that PAH absorption/C is indistinguishable
c from graphite absorption per C (1/3 - 2/3 approx)

            IDIEL_CPE=13
            IDIEL_CPA=30

            CALL INDEX(A,TGR,WAVE,IDIEL_CPA,EPS11,EPS2,ENRE1,ENIM)
            ENRE=1.+ENRE1
            DCXREF=DBLE(ENRE)+(0.D0,1.D0)*DBLE(ENIM)
            XX=DBLE(X)
            CALL MIEV0(XX,DCXREF,PERFCT,MIMCUT,ANYANG,
     &                 NANG,XMU,NMOM,IPOLZN,MOMDIM,PRNT,
     &                 QEXTD,QSCAD,GQSC,PMOM,DCXSFORW,DCXSBACK,S1DCX,
     &                 S2DCX,TFORW,TBACK,SPIKE)
            QABS=(FPAH/3.)*REAL(QEXTD-QSCAD)
            QSCA=(FPAH/3.)*REAL(QSCAD)
            QEXT=(FPAH/3.)*REAL(QEXTD)
            G=(FPAH/3.)*REAL(GQSC)
            G2=(FPAH/3.)*REAL(QSCAD*((PMOM(2,1)+.5D0*PMOM(0,1)))/
     &                       (1.5D0*PMOM(0,1)))
            CALL INDEX(A,TGR,WAVE,IDIEL_CPE,EPS11,EPS2,ENRE1,ENIM)
            ENRE=1.+ENRE1
            DCXREF=DBLE(ENRE)+(0.D0,1.D0)*DBLE(ENIM)
            XX=DBLE(X)
            CALL MIEV0(XX,DCXREF,PERFCT,MIMCUT,ANYANG,
     &                 NANG,XMU,NMOM,IPOLZN,MOMDIM,PRNT,
     &                 QEXTD,QSCAD,GQSC,PMOM,DCXSFORW,DCXSBACK,S1DCX,
     &                 S2DCX,TFORW,TBACK,SPIKE)
            QABS=QABS+(2.*FPAH/3.)*REAL(QEXTD-QSCAD)
            QSCA=QSCA+(2.*FPAH/3.)*REAL(QSCAD)
            QEXT=QEXT+(2.*FPAH/3.)*REAL(QEXTD)
            G=G+(2.*FPAH/3.)*REAL(GQSC)
            G2=G2+(2.*FPAH/3.)*
     &             REAL(QSCAD*((PMOM(2,1)+.5D0*PMOM(0,1)))/
     &                 (1.5D0*PMOM(0,1)))
         ENDIF ! end if wave > .057971 else
      ENDIF ! end if fpah > 0

      IF(NC.GT.NCMAX)THEN
         WRITE(0,*)'Fatal error in qcomp: NC=',NC,' > NCMAX=',NCMAX
         STOP
      ENDIF

      IF(NC.GT.0)THEN
         DO J=1,NC
c*** diagnostic
c            write(0,*)'qcomp_v11 ckpt 80'
c            write(0,*)'about to call index with IDIEL(J)=',IDIEL(J)
c***
            CALL INDEX(A,TGR,WAVE,IDIEL(J),EPS11,EPS2,ENRE1,ENIM)
c*** diagnostic
c            write(0,*)'qcomp_v11 ckpt 90 returned from index'
c            write(0,*)' eps11,eps2=',eps11,eps2
c***
            ENRE=1.+ENRE1
            DCXEPS=1.+EPS11+DCXI*EPS2
            CXINDX=ENRE+DCXI*ENIM

c*** diagnostic
c            write(0,*)'qcomp_v11 ckpt 100, x=',x,' xdipl=',xdipl
c***

            IF(MAGNETIC.AND.X.LT.100.D0)THEN
               RADIUS=DBLE(A)
               LAMBDA=DBLE(WAVE)
               NANGL=181
c*** diagnostic
c               write(0,*)'qcomp_v11 ckpt 110'
c               write(0,*)' about to call magmie'
c               write(0,*)' nangl=',nangl
c               dcxeps=(1.001,0.)
c               write(0,*)' cxeps=',dcxeps
c               write(0,*)' cxmu=',dcxmu
c               write(0,*)' radius=',radius
c               write(0,*)' lambda=',lambda
c***
               CALL MAGMIE(NANGL,DCXEPS,DCXMU,RADIUS,LAMBDA,
     &                     QABSD,QEXTD,QSCAD,GD,
     &                     THETA,S11,S12,S33,S34,POL)
c*** diagnostic
c               write(0,*)'qcomp_v11 ckpt 111'
c               write(0,*)' returned from magmie with'
c               WCABSMAGV=0.75*LAMBDA*QABSD/RADIUS
c               write(0,9100)wcv
c               write(0,*)' qextd=',qextd,' for j=',j
c               write(0,*)' qscad=',qscad
c               write(0,*)' <cos>=',gd
c***
               QQABS(J)=REAL(QABSD)
               QQEXT(J)=REAL(QEXTD)
               QQSCA(J)=REAL(QSCAD)
               GG(J)=REAL(GD)

c evaluate <cos^2>
               G2D=0.D0
               SUM=0.D0
               DO JJ=1,NANGL
                  WGT=DBLE(4-2*MOD(JJ,2))
                  IF(JJ.EQ.1.OR.JJ.EQ.NANGL)WGT=1.D0
                  SN=SIN(PI*THETA(JJ)/180.)
                  CS=COS(PI*THETA(JJ)/180.)
                  G2D=G2D+WGT*S11(JJ)*SN*CS**2
                  SUM=SUM+WGT*S11(JJ)*SN
               ENDDO
               GG2(J)=REAL(G2D/SUM)
c*** diagnostic
c               write(0,*)'qcomp_v11 ckpt 115'
c               write(0,*)' <cos^2>=',gg2(j)
c***
            ELSE

               IF(X.LT.XDIPL)THEN
c*** diagnostic
c                  write(0,*)'qcomp_v11 ckpt 116'
c                  write(0,*)' x=',x
c                  write(0,*)' cxmu=',cxmu
c                  write(0,*)' about to call EMDIPL'
c***
c call emdipl_v2:

                  CALL EMDIPL(X,CXINDX,CXMU,QQABS(J),QQSCA(J),QQEXT(J))

c*** diagnostic
c                  write(0,*)'qcomp_v11 ckpt 117'
c                  write(0,*)' cxmu=',cxmu
c                  write(0,*)' returned from EMDIPL'
c***
                  GG(J)=0.
                  GG2(J)=0.4
               ELSE
                  DCXREF=DBLE(ENRE)+(0.D0,1.D0)*DBLE(ENIM)
                  XX=DBLE(X)
                  Z=X*ABS(CXINDX-1.)
                  Y=X+Z
                  IF(Y.LT.XYMT)THEN

c use Mie theory if x + |m-1|*x < XYMT

c MIEV0 = Wiscombe's routine:
c*** diagnostic
c                     write(0,*)'qcomp_v11 ckpt 120'
c                     write(0,*)' call miev0 for a,wave=',a,wave
c***
                     CALL MIEV0(XX,DCXREF,PERFCT,MIMCUT,ANYANG,
     &                          NANG,XMU,NMOM,IPOLZN,MOMDIM,PRNT,
     &                          QEXTD,QSCAD,GQSC,PMOM,DCXSFORW,DCXSBACK,
     &                          S1DCX,S2DCX,TFORW,TBACK,SPIKE)

c*** diagnostic
c                     write(0,*)'qcomp_v11 ckpt 121'
c                     write(0,*)' returned from miev0'
c***
                     GG(J)=REAL(GQSC/QSCAD)

c PMOM(0,1)= moment of P_0 = 1
c      1,1             P_1 = x
c      2,1             P_2 = 0.5*(3x^2-1)

                     GG2(J)=REAL((PMOM(2,1)+.5D0*PMOM(0,1))/
     &                           (1.5D0*PMOM(0,1)))

c*** diagnostic -- compare Wiscombe's code with bhmie
c                     write(0,*)'qcomp_v11 ckpt 130'
c                     write(0,*)' about to call bhmie'
c                     write(0,*)' with cxindx=',cxindx
c***
c x = real
c cxindx = double complex
c 
                     CALL BHMIE(X,CXINDX,NANG,S1,S2,QEXTPA,QSCAPA,
     &                       QBACK,GPA)
c*** diagnostic
c                     write(0,*)'qcomp_v11 ckpt 131'
c                     write(0,*)' returned from bhmie'
c***
                     if(abs(qextpa/qextd-1.d0).gt.1.d-3)
     &                  write(0,7100)qextd,qextpa,a,wave,cxindx
                     if(abs(qscapa/qscad-1.d0).gt.1.d-3)
     &                  write(0,7200)qscad,qscapa,a,wave,cxindx
c***

                     QQABS(J)=REAL(QEXTD-QSCAD)
                     QQEXT(J)=REAL(QEXTD)
                     QQSCA(J)=REAL(QSCAD)

                  ELSE

                     IF(ABS(DCXREF-1.).LT.M1ADT)THEN

c use anomalous diffraction theory:

                        CALL ADT(XX,DCXREF,QABSD,QSCAD,QEXTD,
     &                           DQDO,THETAS,MXANG,NTHETAS)
                        QQABS(J)=REAL(QABSD)
                        QQEXT(J)=REAL(QEXTD)
                        QQSCA(J)=REAL(QSCAD)
                        GG(J)=1.
                        GG2(J)=1.
                     ELSE

c use geometric optics:

                        CALL GEOMOPT(CXINDX,A,WAVE,QQABS(J),QQSCA(J),
     &                               GG(J))
                        QQEXT(J)=QQABS(J)+QQSCA(J)
                        GG2(J)=0.5+2.*(GG(J)-0.5)**2

                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            G=G+FC(J)*QQSCA(J)*GG(J)
            G2=G2+FC(J)*QQSCA(J)*GG2(J)
            QABS=QABS+FC(J)*QQABS(J)
            QEXT=QEXT+FC(J)*QQEXT(J)
            QSCA=QSCA+FC(J)*QQSCA(J)

c*** diagnostic
c                  write(0,*)'j=',j,' fc(j)=',fc(j),' qqabs(j)=',qqabs(j)
c                  write(0,*)'    qabs=',qabs
c***
         ENDDO
         IF(QSCA.GT.0.)THEN
            G=G/QSCA
            G2=G2/QSCA
         ENDIF
      ENDIF

c 11.12.13 (BTD) now add magnetic contribution to absorption
c 12.12.08 (BTD) modify to suppress this term when a/wave exceeds
c                0.001 (e.g., wave=1cm,   a=10um
c                                  1mm,   a=1um
c                                  100um, a=0.1um
c                (a/wave) -> (a/wave)*[1+(a/wave)/0.001]

      QABS=QABS+(4.D0/3.D0)*(A/WAVE)*WCABSMAGV/(1.+1.D3*(A/WAVE))

      RETURN
 7100 format('Qext=',1PE12.5,',',1PE12.5,' for MIEV0,BHMIE, a=',
     &        1PE12.5,' wave=',1PE12.5,' m=',1P2E10.3)
 7200 format('Qsca=',1PE12.5,',',1PE12.5,' for MIEV0,BHMIE, a=',
     &        1PE12.5,' wave=',1PE12.5,' m=',1P2E10.3)
 9100 format('magnetic wave*Cabs/V=',1pe10.3)
       END
      SUBROUTINE ADT(X,CXREFIN,QABS,QSCA,QEXT,DQDO,
     &               THETAS,MXANG,NTHETAS)
      IMPLICIT NONE

c arguments:

      INTEGER MXANG,NTHETAS
      DOUBLE PRECISION QABS,QEXT,QSCA,X
      DOUBLE PRECISION 
     &   DQDO(0:MXANG),
     &   THETAS(0:MXANG)
      DOUBLE COMPLEX CXREFIN

c local variables:

      INTEGER NLEVMX
      PARAMETER(NLEVMX=28)
      
      INTEGER JLEV,JT,JTH,JU,NLEV,NTERM
      DOUBLE PRECISION BETA,DU,ERT,FAC,FAC2,FAC3,PI,
     &   RHO_0,RHO_1,RHO_2,U,U0,XTH
      DOUBLE COMPLEX CXI,CXIRHO,CXRHO
      DOUBLE COMPLEX
     &   Q(1:NLEVMX),
     &   CXSUM(0:NLEVMX)

c intrinsic functions

      INTRINSIC DIMAG

c external functions

      DOUBLE PRECISION BESJ0
      EXTERNAL BESJ0

c=======================================================================
c Subroutine ADT
c Given:
c   X       = 2*pi*a/lambda
c   CXREFIN = m = complex refractive index
c                 [we adopt convention Im(CXREFIN) > 0 ]
c   THETAS  = scattering angles, JTH=0-NTHETAS
c             NOTE: scattering angles do not need to be uniformly
c                   spaced (or monotonic).
c
c Returns:
c   QABS    = absorption cross section/pi*a^2
c   QSCA    = scattering cross section/pi*a^2
c   QEXT    = extinction cross section/pi*a^2
c   DQDO    = dQ_sca/dOmega for angles JTH=0-NTHETAS
c
c Requires:
c   double precision external function BESJ0(X) to evaluate
c   ordinary Bessel function J_0(x)
c
c Function BESJ0 from netlib SPECFUN library can be used.
c
c Note: parameter NLEVMX limit application of ADT to
c extremely large values of X.  Maximum number of step halvings is
c NLEVMX , so maximum number of integration points is 2**(NLEVMX+1)
c
c-----------------------------------------------------------------------
c
c Subroutine ADT uses anomalous diffraction theory, valid 
c when |m-1| << 1 and  x >> 1 
c
c Note that |m-1|*x can have an arbitrary value
c
c See Chapter 11 of van de Hulst (1957),
c Light Scattering by Small Particles
c
c B.T. Draine, Princeton University Observatory, 2002.05.07
c history
c 02.05.07 (BTD) first written
c 02.05.09 (BTD) calculation of differential scattering cross section
c                totally rewritten to increase computational
c                efficiency by storing and reusing Bessel function
c                values
c 02.07.25 (BTD) concluded that there was no real gain in computational
c                efficiency.  Abandoned strategy to reuse Bessel
c                values.  Scattering angles no longer need to be 
c                uniformly spaced.
c 02.09.06 (BTD) changed NLEVMX from 18 to 20
c 02.12.18 (BTD) changed NLEVMX from 20 to 24
c                changed ERT from 1.D-3 to 1.D-4
c                changed NLEVMX from 24 to 26
c end history
c=======================================================================
      DATA CXI/(0.D0,1.D0)/,ERT/1.D-4/
      SAVE CXI,ERT
c
c Following four statements should be enabled if NOT using g77.
c They assume that the compiler supports double complex, since the
c statements DBLE and DIMAG are used.  If double complex is not available
c (see above) you will need to change DIMAG to AIMAG
c
c If using g77, following statments could be commented out, as 
c REALPART and IMAGPART are g77 intrinsic functions
c However, they do not need to be commented out.

      DOUBLE PRECISION IMAGPART,REALPART
      DOUBLE COMPLEX DCXVAR
      REALPART(DCXVAR)=(DBLE(DCXVAR))
      IMAGPART(DCXVAR)=(DIMAG(DCXVAR))
c=======================================================================      

      PI=4.D0*ATAN(1.D0)
      CXRHO=2.D0*X*(CXREFIN-1.D0)
      CXIRHO=CXI*CXRHO
      RHO_0=ABS(CXRHO)
      RHO_1=REALPART(CXRHO)
      RHO_2=IMAGPART(CXRHO)
      IF(ABS(RHO_1).GT.0.D0)THEN
         BETA=ATAN(RHO_2/RHO_1)
      ELSE
         IF(RHO_2.GT.0.D0)THEN
            BETA=0.5D0*PI
         ELSE
            BETA=-0.5D0*PI
         ENDIF
      ENDIF
      IF(RHO_0.LT.1.D-3)THEN
         QEXT=(4.D0/3.D0)*RHO_2+0.5D0*(RHO_1**2-RHO_2**2)
         QABS=(4.D0/3.D0)*RHO_2-RHO_2**2
         QSCA=0.5D0*RHO_0**2
      ELSE
         FAC=EXP(-RHO_2)
         FAC2=FAC*FAC
         QEXT=2.D0+4.D0*(COS(2.D0*BETA)-
     &        FAC*(COS(RHO_1-2.D0*BETA)-RHO_0*SIN(RHO_1-BETA)))/RHO_0**2
         QABS=1.D0+(FAC2+0.5D0*(FAC2-1.D0)/RHO_2)/RHO_2
         QSCA=QEXT-QABS
      ENDIF

c-----------------------------------------------------------------------
c calculate differential scattering cross section
c
      FAC=X*X/PI

      DO JTH=0,NTHETAS
         XTH=X*THETAS(JTH)

         IF(XTH.EQ.0.D0)THEN

c if scattering angle is zero, evaluate forward scattering analytically:

            DQDO(JTH)=FAC*
     &                ABS(0.5D0+(1.D0+(CXIRHO-1.D0)*EXP(CXIRHO))/
     &                         CXRHO**2)**2
         ELSE

c evaluate using adaptive Simpson's rule, using U=cos(tau) as integration
c variable.  Note that integrand vanishes at endpoints U=0 and U=1.

            U=0.5D0
            CXSUM(0)=(1.D0-EXP(CXIRHO*SQRT(1.D0-U*U)))*BESJ0(XTH*U)*U
            DO JLEV=1,NLEVMX
               NTERM=2**JLEV
               DU=1.D0/DBLE(NTERM)
               U0=-0.5D0*DU
               CXSUM(JLEV)=0.D0
               DO JU=1,NTERM
                  U=U0+DU*JU
                  FAC2=BESJ0(XTH*U)*U
                  FAC3=SQRT(1.D0-U*U)
                  CXSUM(JLEV)=CXSUM(JLEV)+(1.D0-EXP(CXIRHO*FAC3))*FAC2
               ENDDO

c quadrature using Simpson's rule:

               Q(JLEV)=2.D0*CXSUM(JLEV)
               DO JT=0,JLEV-1
                  Q(JLEV)=Q(JLEV)+CXSUM(JT)
               ENDDO
               Q(JLEV)=Q(JLEV)*DU/3.D0
               NLEV=JLEV

c check to see whether Simpson's rule integration has converged

               IF(ABS(Q(NLEV)-Q(NLEV-1))/ABS(Q(NLEV)+Q(NLEV-1)).LT.ERT)
     &            GOTO 3500
            ENDDO

c if reach here, we failed to meet error tolerance.
c write out warning, but use best estimate for diff. scatt. cross section

            WRITE(0,9200)NLEVMX
 3500       DQDO(JTH)=FAC*ABS(Q(NLEV))**2
         ENDIF
      ENDDO
      
      RETURN

 9200 FORMAT('Warning: adt failed to meet error tolerance in NLEV=',I2,
     &          ' doublings')
      END




      SUBROUTINE ALPHA_EDDY(IDIEL,OMEGA,RADIUS,TGR,CXALPHAV_EDDY)
      IMPLICIT NONE
c arguments
      INTEGER IDIEL
      DOUBLE PRECISION OMEGA,RADIUS,TGR
      DOUBLE COMPLEX CXALPHAV_EDDY

c local variables

      REAL ENIM_SP,ENRE1_SP,EPS11_SP,EPS2_SP,RADIUS_SP,
     &   TGR_SP,WAVE_SP
      DOUBLE PRECISION PI
      COMPLEX CXI,CXY
      COMPLEX SIN,COS
c-----------------------------------------------------------------------
c subroutine ALPHA_EDDY
c given:
c   IDIEL  = composition identifier
c   OMEGA  = frequency (rad/s)
c   RADIUS = radius (um)
c returns
c   CXALPHAV_EDDY = alpha/Volume (dimensionless)
c            where alpha = magnetic polarizability due to eddy currents
c            computed following Landau & Lifshitz
c B.T. Draine, Princeton Univ. Observatory, 2011.12.08
c history
c 2011.12.08 (BTD) first written
c end history
c-----------------------------------------------------------------------
      CXI=(0.D0,1.D0)
      PI=4.D0*ATAN(1.D0)
      RADIUS_SP=REAL(RADIUS)
      TGR_SP=REAL(TGR)
      WAVE_SP=2.*PI*2.997925E14/OMEGA
      CALL INDEX(RADIUS_SP,TGR_SP,WAVE_SP,IDIEL,EPS11_SP,EPS2_SP,
     &           ENRE1_SP,ENIM_SP)
      CXY=(OMEGA*RADIUS/2.997925D14)*(1.+ENRE1_SP+CXI*ENIM_SP)
c*** diagnostic
c      write(7,fmt='(a,1pe10.3,a,1pe10.3,a,1pe10.3,a,1pe10.2,1pe10.3)')
c     &        'a=',radius,' n=',(1.+enre1_sp),' k=',enim_sp,' y=',cxy
c***
      IF(ABS(CXY).GT.0.1)THEN
c*** diagnostic
c         WRITE(7,fmt='(A,1pe10.3)')'eddy calc, |y|=',abs(cxy)
c***
         CXALPHAV_EDDY=(3./(8.*PI))*
     &                 (3./CXY**2-3.*COS(CXY)/(CXY*SIN(CXY))-1.)
      ELSE
         CXALPHAV_EDDY=CXY**2*(1.+(2./21.)*CXY**2+CXY**4/105.)/(40.*PI)
      ENDIF
      RETURN
      END
      subroutine amuelr (s1,s2,theta,numang,sign, 
     &                                     s11nor,s11,s12,s33,s34,pol)
c
c **********************************************
c
c subroutine amuelr computes Mueller matrix elements for either
c spheres or infinite cylinders
c
c       Merrill Milham      >>> version: 2.0 <<<        JANUARY 1994
c
c inputs :
c
c           s1 = amplitude scattering matrix element array.(complex*16)
c           s2 = amplitude scattering matrix element array.(complex*16)
c        theta = scattering angle array (real*8)
c       numang = the number of scattering angles, i.e., number of
c                elements in s1, s2. (integer)
c         sign = arbitrary sign (+1 or -1) used to adjust the sign of
c                s34 according to the users convention. (integer)
c
c outputs :
c
c      s11nor = s11 for a scattering angle of zero degrees, which
c               is used to normalize s11(theta) (real*8)
c         s11 = Mueller matrix element 1,l (real*8)
c         s12 = Mueller matrix element 1,2 (real*8)
c         s33 = Mueller matrix element 3,3 (real*8)
c         s34 = Mueller matrix element 3,4 (real*8)
c         pol = polarization = pol=-s12/s11 (real*8)
c
c         s12,s33,s34 are normalized by s11(theta)
c
c subroutines used: none
c
c *************************************************
c
      implicit none
*
      complex*16 s1(1),s2(1)
      real*8 theta(1),s11nor
      integer numang,sign
*
      real*8 s11(1),s12(1),s33(1),s34(1)
      real*8 pol(1)
*
      real*8 rs1,is1,rs2,is2,s11i
      real*8 ts11,ts12,half,zero,one
      parameter (half=0.5d0,zero=0.d0,one=1.d0)
      complex*16 ts
      integer i
c
      do 100 i=1,numang
      rs1=dble(s1(i))
      is1=dimag(s1(i))
      rs2=dble(s2(i))
c
c page 37 
c
      is2=dimag(s2(i)) 
*
      ts11=rs1*rs1+is1*is1
      ts12=rs2*rs2+is2*is2
      s11(i)=half*(ts11+ts12)
      s12(i)=half*(ts12-ts11) 
*
      if(s11(i).ne.zero) then
            s11i=s11(i)
                          else
            s11i=one
            write(*,*)
            write(*,*) 'Unnormalized Mueller matrix elements for ',
     &                  theta(i),' deg.'
            write(*, *)
      end if 
*
      pol(i)=-s12(i)/s11i 
*
      ts=s2(i)*dconjg(s1(i))
      s33(i)=dble(ts)
      s33(i)=s33(i)/s11i
      s34(i)=sign*dimag(ts)
      s34(i)=s34(i)/s11i 
*
      if(i.eq.1) s11nor=s11i 
*
      s11(i)=s11i/s11nor
100   continue
c
      return
c
      end
      SUBROUTINE BHMIE(X,REFREL,NANG,S1,S2,QEXT,QSCA,QBACK,GSCA)
      IMPLICIT NONE

C Declare parameters:
C Note: important that MXNANG be consistent with dimension of S1 and S2
C       in calling routine!

      INTEGER MXNANG,NMXX
C      PARAMETER(MXNANG=1000,NMXX=15000)
      PARAMETER(MXNANG=1000,NMXX=150000)

C Arguments:

      INTEGER NANG
      REAL GSCA,QBACK,QEXT,QSCA,X
      COMPLEX REFREL
      COMPLEX S1(2*MXNANG-1),S2(2*MXNANG-1)

C Local variables:

      LOGICAL SINGLE
      INTEGER J,JJ,N,NSTOP,NMX,NN
      DOUBLE PRECISION CHI,CHI0,CHI1,DANG,DX,EN,FN,P,PII,PSI,PSI0,PSI1,
     &                 THETA,XSTOP,YMOD
      DOUBLE PRECISION
     &   AMU(MXNANG),
     &   PI(MXNANG),
     &   PI0(MXNANG),
     &   PI1(MXNANG),
     &   TAU(MXNANG)
      DOUBLE COMPLEX
     &   DCXS1(2*MXNANG-1),
     &   DCXS2(2*MXNANG-1)

      INTRINSIC DCMPLX,DIMAG

C***********************************************************************
C
C Subroutine BHMIE is derived from the Bohren-Huffman Mie scattering
C     subroutine to calculate scattering and absorption by a homogenous
C     isotropic sphere.
C Given:
C    X = 2*pi*a/lambda
C    REFREL = (complex refr. index of sphere)/(real index of medium)
C    NANG = number of angles between 0 and 90 degrees
C           (will calculate 2*NANG-1 directions from 0 to 180 deg.)
C           if called with NANG<2, will set NANG=2 and will compute
C           scattering for theta=0,90,180.
C Returns:
C    S1(1 - 2*NANG-1) = -i*f_22 (incid. E perp. to scatt. plane,
C                                scatt. E perp. to scatt. plane)
C    S2(1 - 2*NANG-1) = -i*f_11 (incid. E parr. to scatt. plane,
C                                scatt. E parr. to scatt. plane)
C    QEXT = C_ext/pi*a**2 = efficiency factor for extinction
C    QSCA = C_sca/pi*a**2 = efficiency factor for scattering
C    QBACK = 4.*pi*(dC_sca/domega)/pi*a**2
C          = backscattering efficiency
C    GSCA = <cos(theta)> for scattering
C
C S1 and S2 are the diagonal elements of the "amplitude scattering matrix"
C (see eq. 3.12 of Bohren & Huffman 1983) -- the off-diagonal elements
C vanish for a spherical target.
C For unpolarized incident light, the intensity of scattered light a
C distance r from the sphere is just
C          1
C  I_s = ------ * I_in * S_11
C        (kr)^2
C
C where k=2*pi/lambda 
C and the "Muller matrix element" S_11 = 0.5*( |S_1|^2 + |S_2|^2 )
C
C for incident light polarized perp to the scattering plane,
C the scattered light is polarized perp to the scattering plane
C with intensity I_s = I_in * |S_1|^2 / (kr)^2
C
C for incident light polarized parallel to the scattering plane,
C the scattered light is polarized parallel to the scattering plane
C with intensity I_s = I_in * |S_2|^2 / (kr)^2
C
C History:
C Original program taken from Bohren and Huffman (1983), Appendix A
C Modified by B.T.Draine, Princeton Univ. Obs., 90.10.26
C in order to compute <cos(theta)>
C 91.05.07 (BTD): Modified to allow NANG=1
C 91.08.15 (BTD): Corrected error (failure to initialize P)
C 91.08.15 (BTD): Modified to enhance vectorizability.
C 91.08.15 (BTD): Modified to make NANG=2 if called with NANG=1
C 91.08.15 (BTD): Changed definition of QBACK.
C 92.01.08 (BTD): Converted to full double precision and double complex
C                 eliminated 2 unneed lines of code
C                 eliminated redundant variables (e.g. APSI,APSI0)
C                 renamed RN -> EN = double precision N
C                 Note that DOUBLE COMPLEX and DCMPLX are not part
C                 of f77 standard, so this version may not be fully
C                 portable.  In event that portable version is
C                 needed, use src/bhmie_f77.f
C 93.06.01 (BTD): Changed AMAX1 to generic function MAX
C 98.09.17 (BTD): Added variable "SINGLE" and warning in event that
C                 code is used with single-precision arithmetic (i.e.,
C                 compiler does not support DOUBLE COMPLEX)
C 99.02.17 (BTD): Replaced calls to REAL() and IMAG() by
C                 REALPART() and IMAGPART() for compatibility with g77
C                 Note that when code is used with standard f77 
C                 compilers, it is now necessary to enable two lines
C                 defining functions REALPART(X) and IMAGPART(X)
C 99.02.19 (BTD): added lines to be enabled to properly define
C                 REALPART() and IMAGPART() if NOT using g77
C                 ***see below!!***
C 01.02.16 (BTD): added IMPLICIT NONE
C 01.02.27 (BTD): changed definition of QBACK back to convention of
C                 Bohren & Huffman and others:
C                 Q_back = 4.*pi*(dC_sca/dOmega)/(pi*a^2) in backward
C                          direction
c 02.03.09 (BTD): defined statement function REALPART_SP to
c                 avoid warning regarding type conversion when taking
c                 real part of S1(1) to evaluate QEXT
c                 some cleanup regarding type conversion
c 02.05.30 (BTD): introduced internal double complex arrays DCXS1,DCXS2
c                 to possibly increase accuracy during summations.
c                 After summations, output scattering amplitudes
c                 via single complex arrays S1,S2 as before.
c                 Usage of this routine is unaffected by change.
c                 Note: no longer need statement function REALPART_SP
c 02.09.18 (BTD): Error in evaluation of QBACK reported by Okada Yasuhiko
c                 Was calculating QBACK using S1 rather than DCXS1
c                 Corrected.
c 02.10.16 (BTD): Added comments explaining definition of S_1 and S_2 .
C end history
C
C***********************************************************************
C 
C This module is dependent on whether compiler supports double precision
C complex variables:
C
C If your compiler does NOT support double complex, comment out following
C three lines, and uncomment corresponding 3 lines further below
C
      DOUBLE COMPLEX AN,AN1,BN,BN1,DREFRL,XI,XI1,Y
      DOUBLE COMPLEX D(NMXX)
      PARAMETER(SINGLE=.FALSE.)

C      COMPLEX AN,AN1,BN,BN1,DREFRL,XI,XI1,Y
C      COMPLEX D(NMXX)
C      PARAMETER(SINGLE=.TRUE.)

C**********************************************************************

C Following five statements should be enabled if NOT using g77.
C They assume that the compiler supports double complex, since the
C statements DBLE and DIMAG are used.  If double complex is not available
C (see above) you will need to change DIMAG to AIMAG
C
C If using g77, following statements could be commented out, as 
C REALPART and IMAGPART are g77 intrinsic functions
C However, they do not need to be commented out.

      DOUBLE COMPLEX DPCX
      DOUBLE PRECISION REALPART
      DOUBLE PRECISION IMAGPART
      REALPART(DPCX)=(DBLE(DPCX))
      IMAGPART(DPCX)=(DIMAG(DPCX))
      
C***********************************************************************
c*** diagnostic
c      write(0,*)'bhmie ckpt 1'
c      write(0,*)' refrel=',refrel
c***
C*** Safety checks

      IF(SINGLE)WRITE(0,*)'Warning: this version of bhmie uses only ',
     &          'single precision complex numbers!'
      IF(NANG.GT.MXNANG)STOP'***Error: NANG > MXNANG in bhmie'
      IF(NANG.LT.2)NANG=2

C*** Obtain pi:

      PII=4.D0*ATAN(1.D0)
      DX=X
      DREFRL=REFREL
      Y=X*DREFRL
      YMOD=ABS(Y)

C*** Series expansion terminated after NSTOP terms
C    Logarithmic derivatives calculated from NMX on down
      XSTOP=X+4.*X**0.3333+2.
      NMX=NINT(MAX(XSTOP,YMOD))+15
c*** diagnostic
c      write(0,*)'bhmie ckpt 50'
c      write(0,*)' x=',x
c      write(0,*)' refrel=',refrel
c      write(0,*)' drefrl=',drefrl
c      write(0,*)' y=',y
c      write(0,*)' xstop=',xstop
c      write(0,*)' ymod=',ymod
c      write(0,*)' nmx=',nmx
c***
C BTD experiment 91.1.15: add one more term to series and compare results
C      NMX=MAX(XSTOP,YMOD)+16
C test: compute 7001 wavelengths between .0001 and 1000 micron
C for a=1.0micron SiC grain.  When NMX increased by 1, only a single
C computed number changed (out of 4*7001) and it only changed by 1/8387
C conclusion: we are indeed retaining enough terms in series!

      NSTOP=NINT(XSTOP)

      IF(NMX.GT.NMXX)THEN
         WRITE(0,*)'Error: NMX > NMXX=',NMXX,' for |m|x=',YMOD
         STOP
      ENDIF

C*** Require NANG.GE.1 in order to calculate scattering intensities

      DANG=0.
      IF(NANG.GT.1)DANG=.5*PII/DBLE(NANG-1)
      DO J=1,NANG
         THETA=DBLE(J-1)*DANG
         AMU(J)=COS(THETA)
      ENDDO
      DO J=1,NANG
         PI0(J)=0.
         PI1(J)=1.
      ENDDO
      NN=2*NANG-1
      DO J=1,NN
         DCXS1(J)=(0.D0,0.D0)
         DCXS2(J)=(0.D0,0.D0)
      ENDDO

C*** Logarithmic derivative D(J) calculated by downward recurrence
C    beginning with initial value (0.,0.) at J=NMX

c***  diagnostic
c      write(0,*)'bhmie ckpt 100'
c      write(0,*)' nn=',nn,' nmx=',nmx
c***
      D(NMX)=(0.,0.)
      NN=NMX-1

c*** diagnostic
c      write(0,*)'bhmie ckpt 101'
c***
      DO N=1,NN
         EN=NMX-N+1
         D(NMX-N)=(EN/Y)-(1./(D(NMX-N+1)+EN/Y))
      ENDDO
c*** diagnostic
c      write(0,*)'bhmie ckpt 110'
c***
C*** Riccati-Bessel functions with real argument X
C    calculated by upward recurrence

      PSI0=COS(DX)
      PSI1=SIN(DX)
      CHI0=-SIN(DX)
      CHI1=COS(DX)
      XI1=DCMPLX(PSI1,-CHI1)
      QSCA=0.E0
      GSCA=0.E0
      P=-1.
      DO N=1,NSTOP
         EN=N
         FN=(2.E0*EN+1.)/(EN*(EN+1.))

C for given N, PSI  = psi_n        CHI  = chi_n
C              PSI1 = psi_{n-1}    CHI1 = chi_{n-1}
C              PSI0 = psi_{n-2}    CHI0 = chi_{n-2}
C Calculate psi_n and chi_n

         PSI=(2.E0*EN-1.)*PSI1/DX-PSI0
         CHI=(2.E0*EN-1.)*CHI1/DX-CHI0
         XI=DCMPLX(PSI,-CHI)

C*** Store previous values of AN and BN for use
C    in computation of g=<cos(theta)>

         IF(N.GT.1)THEN
            AN1=AN
            BN1=BN
         ENDIF

C*** Compute AN and BN:

         AN=(D(N)/DREFRL+EN/DX)*PSI-PSI1
         AN=AN/((D(N)/DREFRL+EN/DX)*XI-XI1)
         BN=(DREFRL*D(N)+EN/DX)*PSI-PSI1
         BN=BN/((DREFRL*D(N)+EN/DX)*XI-XI1)

C*** Augment sums for Qsca and g=<cos(theta)>

         QSCA=QSCA+REAL((2.*EN+1.)*(ABS(AN)**2+ABS(BN)**2))
         GSCA=GSCA+REAL(((2.*EN+1.)/(EN*(EN+1.)))*
     &        (REALPART(AN)*REALPART(BN)+IMAGPART(AN)*IMAGPART(BN)))
         IF(N.GT.1)THEN
            GSCA=GSCA+REAL(((EN-1.)*(EN+1.)/EN)*
     &      (REALPART(AN1)*REALPART(AN)+IMAGPART(AN1)*IMAGPART(AN)+
     &      REALPART(BN1)*REALPART(BN)+IMAGPART(BN1)*IMAGPART(BN)))
         ENDIF

C*** Now calculate scattering intensity pattern
C    First do angles from 0 to 90

         DO J=1,NANG
            JJ=2*NANG-J
            PI(J)=PI1(J)
            TAU(J)=EN*AMU(J)*PI(J)-(EN+1.)*PI0(J)
            DCXS1(J)=DCXS1(J)+FN*(AN*PI(J)+BN*TAU(J))
            DCXS2(J)=DCXS2(J)+FN*(AN*TAU(J)+BN*PI(J))
         ENDDO

C*** Now do angles greater than 90 using PI and TAU from
C    angles less than 90.
C    P=1 for N=1,3,...; P=-1 for N=2,4,...

         P=-P
         DO J=1,NANG-1
            JJ=2*NANG-J
            DCXS1(JJ)=DCXS1(JJ)+FN*P*(AN*PI(J)-BN*TAU(J))
            DCXS2(JJ)=DCXS2(JJ)+FN*P*(BN*PI(J)-AN*TAU(J))
         ENDDO
         PSI0=PSI1
         PSI1=PSI
         CHI0=CHI1
         CHI1=CHI
         XI1=DCMPLX(PSI1,-CHI1)

C*** Compute pi_n for next value of n
C    For each angle J, compute pi_n+1
C    from PI = pi_n , PI0 = pi_n-1

         DO J=1,NANG
            PI1(J)=((2.*EN+1.)*AMU(J)*PI(J)-(EN+1.)*PI0(J))/EN
            PI0(J)=PI(J)
         ENDDO
      ENDDO

C*** Have summed sufficient terms.
C    Now compute QSCA,QEXT,QBACK,and GSCA

      GSCA=REAL(2.D0*GSCA/QSCA)
      QSCA=REAL((2.D0/(DX*DX))*QSCA)
      QEXT=REAL((4.D0/(DX*DX))*REALPART(DCXS1(1)))
      QBACK=REAL(4.D0*(ABS(DCXS1(2*NANG-1))/DX)**2)

C prepare single precision complex scattering amplitude for output

      DO J=1,2*NANG-1
         S1(J)=CMPLX(DCXS1(J))
         S2(J)=CMPLX(DCXS2(J))
      ENDDO

      RETURN
      END
*DECK D1MACH
      DOUBLE PRECISION FUNCTION D1MACH (I)
C***BEGIN PROLOGUE  D1MACH
C***PURPOSE  Return floating point machine dependent constants.
C***LIBRARY   SLATEC
C***CATEGORY  R1
C***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Fox, P. A., (Bell Labs)
C           Hall, A. D., (Bell Labs)
C           Schryer, N. L., (Bell Labs)
C***DESCRIPTION
C
C   D1MACH can be used to obtain machine-dependent parameters for the
C   local machine environment.  It is a function subprogram with one
C   (input) argument, and can be referenced as follows:
C
C        D = D1MACH(I)
C
C   where I=1,...,5.  The (output) value of D above is determined by
C   the (input) value of I.  The results for various values of I are
C   discussed below.
C
C   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
C   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C   D1MACH( 3) = B**(-T), the smallest relative spacing.
C   D1MACH( 4) = B**(1-T), the largest relative spacing.
C   D1MACH( 5) = LOG10(B)
C
C   Assume double precision numbers are represented in the T-digit,
C   base-B form
C
C              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
C   EMIN .LE. E .LE. EMAX.
C
C   The values of B, T, EMIN and EMAX are provided in I1MACH as
C   follows:
C   I1MACH(10) = B, the base.
C   I1MACH(14) = T, the number of base-B digits.
C   I1MACH(15) = EMIN, the smallest exponent E.
C   I1MACH(16) = EMAX, the largest exponent E.
C
C   To alter this function for a particular environment, the desired
C   set of DATA statements should be activated by removing the C from
C   column 1.  Also, the values of D1MACH(1) - D1MACH(4) should be
C   checked for consistency with the local operating system.
C
C***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
C                 a portable library, ACM Transactions on Mathematical
C                 Software 4, 2 (June 1978), pp. 177-188.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   750101  DATE WRITTEN
C   890213  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900618  Added DEC RISC constants.  (WRB)
C   900723  Added IBM RS 6000 constants.  (WRB)
C   900911  Added SUN 386i constants.  (WRB)
C   910710  Added HP 730 constants.  (SMR)
C   911114  Added Convex IEEE constants.  (WRB)
C   920121  Added SUN -r8 compiler option constants.  (WRB)
C   920229  Added Touchstone Delta i860 constants.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   920625  Added CONVEX -p8 and -pd8 compiler option constants.
C           (BKS, WRB)
C   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
C***END PROLOGUE  D1MACH
C
ccc      INTEGER SMALL(4)
ccc      INTEGER LARGE(4)
ccc      INTEGER RIGHT(4)
ccc      INTEGER DIVER(4)
ccc      INTEGER LOG10(4)

      DOUBLE PRECISION SMALL(4)
      DOUBLE PRECISION LARGE(4)
      DOUBLE PRECISION RIGHT(4)
      DOUBLE PRECISION DIVER(4)
      DOUBLE PRECISION LOG10(4)
      
C
      DOUBLE PRECISION DMACH(5)
      SAVE DMACH
C
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FDFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA SMALL(1), SMALL(2) / 16#00100000, 16#00000000 /
C     DATA LARGE(1), LARGE(2) / 16#7FFFFFFF, 16#FFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / 16#3CA00000, 16#00000000 /
C     DATA DIVER(1), DIVER(2) / 16#3CB00000, 16#00000000 /
C     DATA LOG10(1), LOG10(2) / 16#3FD34413, 16#509F79FF /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA SMALL(1) / ZC00800000 /
C     DATA SMALL(2) / Z000000000 /
C     DATA LARGE(1) / ZDFFFFFFFF /
C     DATA LARGE(2) / ZFFFFFFFFF /
C     DATA RIGHT(1) / ZCC5800000 /
C     DATA RIGHT(2) / Z000000000 /
C     DATA DIVER(1) / ZCC6800000 /
C     DATA DIVER(2) / Z000000000 /
C     DATA LOG10(1) / ZD00E730E7 /
C     DATA LOG10(2) / ZC77800DC0 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
C
C     DATA SMALL(1) / O1771000000000000 /
C     DATA SMALL(2) / O0000000000000000 /
C     DATA LARGE(1) / O0777777777777777 /
C     DATA LARGE(2) / O0007777777777777 /
C     DATA RIGHT(1) / O1461000000000000 /
C     DATA RIGHT(2) / O0000000000000000 /
C     DATA DIVER(1) / O1451000000000000 /
C     DATA DIVER(2) / O0000000000000000 /
C     DATA LOG10(1) / O1157163034761674 /
C     DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
C
C     DATA SMALL(1) / O1771000000000000 /
C     DATA SMALL(2) / O7770000000000000 /
C     DATA LARGE(1) / O0777777777777777 /
C     DATA LARGE(2) / O7777777777777777 /
C     DATA RIGHT(1) / O1461000000000000 /
C     DATA RIGHT(2) / O0000000000000000 /
C     DATA DIVER(1) / O1451000000000000 /
C     DATA DIVER(2) / O0000000000000000 /
C     DATA LOG10(1) / O1157163034761674 /
C     DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA SMALL(1) / Z"3001800000000000" /
C     DATA SMALL(2) / Z"3001000000000000" /
C     DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" /
C     DATA LARGE(2) / Z"4FFE000000000000" /
C     DATA RIGHT(1) / Z"3FD2800000000000" /
C     DATA RIGHT(2) / Z"3FD2000000000000" /
C     DATA DIVER(1) / Z"3FD3800000000000" /
C     DATA DIVER(2) / Z"3FD3000000000000" /
C     DATA LOG10(1) / Z"3FFF9A209A84FBCF" /
C     DATA LOG10(2) / Z"3FFFF7988F8959AC" /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA SMALL(1) / 00564000000000000000B /
C     DATA SMALL(2) / 00000000000000000000B /
C     DATA LARGE(1) / 37757777777777777777B /
C     DATA LARGE(2) / 37157777777777777777B /
C     DATA RIGHT(1) / 15624000000000000000B /
C     DATA RIGHT(2) / 00000000000000000000B /
C     DATA DIVER(1) / 15634000000000000000B /
C     DATA DIVER(2) / 00000000000000000000B /
C     DATA LOG10(1) / 17164642023241175717B /
C     DATA LOG10(2) / 16367571421742254654B /
C
C     MACHINE CONSTANTS FOR THE CELERITY C1260
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fn OR -pd8 COMPILER OPTION
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CC0000000000000' /
C     DATA DMACH(4) / Z'3CD0000000000000' /
C     DATA DMACH(5) / Z'3FF34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fi COMPILER OPTION
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -p8 COMPILER OPTION
C
C     DATA DMACH(1) / Z'00010000000000000000000000000000' /
C     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3F900000000000000000000000000000' /
C     DATA DMACH(4) / Z'3F910000000000000000000000000000' /
C     DATA DMACH(5) / Z'3FFF34413509F79FEF311F12B35816F9' /
C
C     MACHINE CONSTANTS FOR THE CRAY
C
C     DATA SMALL(1) / 201354000000000000000B /
C     DATA SMALL(2) / 000000000000000000000B /
C     DATA LARGE(1) / 577767777777777777777B /
C     DATA LARGE(2) / 000007777777777777774B /
C     DATA RIGHT(1) / 376434000000000000000B /
C     DATA RIGHT(2) / 000000000000000000000B /
C     DATA DIVER(1) / 376444000000000000000B /
C     DATA DIVER(2) / 000000000000000000000B /
C     DATA LOG10(1) / 377774642023241175717B /
C     DATA LOG10(2) / 000007571421742254654B /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
C     STATIC DMACH(5)
C
C     DATA SMALL /    20K, 3*0 /
C     DATA LARGE / 77777K, 3*177777K /
C     DATA RIGHT / 31420K, 3*0 /
C     DATA DIVER / 32020K, 3*0 /
C     DATA LOG10 / 40423K, 42023K, 50237K, 74776K /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING G_FLOAT
C
C     DATA DMACH(1) / '0000000000000010'X /
C     DATA DMACH(2) / 'FFFFFFFFFFFF7FFF'X /
C     DATA DMACH(3) / '0000000000003CC0'X /
C     DATA DMACH(4) / '0000000000003CD0'X /
C     DATA DMACH(5) / '79FF509F44133FF3'X /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING IEEE_FORMAT
C
C     DATA DMACH(1) / '0010000000000000'X /
C     DATA DMACH(2) / '7FEFFFFFFFFFFFFF'X /
C     DATA DMACH(3) / '3CA0000000000000'X /
C     DATA DMACH(4) / '3CB0000000000000'X /
C     DATA DMACH(5) / '3FD34413509F79FF'X /
C
C     MACHINE CONSTANTS FOR THE DEC RISC
C
C     DATA SMALL(1), SMALL(2) / Z'00000000', Z'00100000'/
C     DATA LARGE(1), LARGE(2) / Z'FFFFFFFF', Z'7FEFFFFF'/
C     DATA RIGHT(1), RIGHT(2) / Z'00000000', Z'3CA00000'/
C     DATA DIVER(1), DIVER(2) / Z'00000000', Z'3CB00000'/
C     DATA LOG10(1), LOG10(2) / Z'509F79FF', Z'3FD34413'/
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     USING D_FLOATING
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
C
C     DATA SMALL(1), SMALL(2) /        128,           0 /
C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
C     DATA RIGHT(1), RIGHT(2) /       9344,           0 /
C     DATA DIVER(1), DIVER(2) /       9472,           0 /
C     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /
C
C     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /
C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     USING G_FLOATING
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
C
C     DATA SMALL(1), SMALL(2) /         16,           0 /
C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
C     DATA RIGHT(1), RIGHT(2) /      15552,           0 /
C     DATA DIVER(1), DIVER(2) /      15568,           0 /
C     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /
C
C     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /
C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C     (ASSUMING REAL*8 IS THE DEFAULT DOUBLE PRECISION)
C
C     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X /
C     DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X /
C     DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X /
C     DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X /
C     DATA LOG10(1), LOG10(2) / '3FD34413'X,'509F79FF'X /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
C     DATA LARGE(1), LARGE(2) / '37777777, '37777577 /
C     DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 /
C     DATA DIVER(1), DIVER(2) / '20000000, '00000334 /
C     DATA LOG10(1), LOG10(2) / '23210115, '10237777 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 /
C     DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 /
C     DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 /
C     DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 /
C     DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 /
C
C     MACHINE CONSTANTS FOR THE HP 730
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     THREE WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 /
C     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B /
C     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B /
C     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B /
C     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) /  40000B,       0 /
C     DATA SMALL(3), SMALL(4) /       0,       1 /
C     DATA LARGE(1), LARGE(2) /  77777B, 177777B /
C     DATA LARGE(3), LARGE(4) / 177777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) /  40000B,       0 /
C     DATA RIGHT(3), RIGHT(4) /       0,    225B /
C     DATA DIVER(1), DIVER(2) /  40000B,       0 /
C     DATA DIVER(3), DIVER(4) /       0,    227B /
C     DATA LOG10(1), LOG10(2) /  46420B,  46502B /
C     DATA LOG10(3), LOG10(4) /  76747B, 176377B /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B /
C     DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B /
C     DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B /
C     DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B /
C     DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 /
C     DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C     ASSUMES THAT ALL ARITHMETIC IS DONE IN DOUBLE PRECISION
C     ON 8088, I.E., NOT IN 80 BIT FORM FOR THE 8087.
C
      DATA SMALL(1) / 2.23D-308  /
      DATA LARGE(1) / 1.79D+308  /
      DATA RIGHT(1) / 1.11D-16   /
      DATA DIVER(1) / 2.22D-16   /
      DATA LOG10(1) / 0.301029995663981195D0 /
C
C     MACHINE CONSTANTS FOR THE IBM RS 6000
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE INTEL i860
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
C
C     DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 /
C     DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 /
C     DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 /
C     DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 /
C     DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
C
C     DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 /
C     DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 /
C     DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 /
C     DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 /
C     DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /    8388608,           0 /
C     DATA LARGE(1), LARGE(2) / 2147483647,          -1 /
C     DATA RIGHT(1), RIGHT(2) /  612368384,           0 /
C     DATA DIVER(1), DIVER(2) /  620756992,           0 /
C     DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 /
C
C     DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 /
C     DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 /
C     DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 /
C     DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 /
C     DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /    128,      0 /
C     DATA SMALL(3), SMALL(4) /      0,      0 /
C     DATA LARGE(1), LARGE(2) /  32767,     -1 /
C     DATA LARGE(3), LARGE(4) /     -1,     -1 /
C     DATA RIGHT(1), RIGHT(2) /   9344,      0 /
C     DATA RIGHT(3), RIGHT(4) /      0,      0 /
C     DATA DIVER(1), DIVER(2) /   9472,      0 /
C     DATA DIVER(3), DIVER(4) /      0,      0 /
C     DATA LOG10(1), LOG10(2) /  16282,   8346 /
C     DATA LOG10(3), LOG10(4) / -31493, -12296 /
C
C     DATA SMALL(1), SMALL(2) / O000200, O000000 /
C     DATA SMALL(3), SMALL(4) / O000000, O000000 /
C     DATA LARGE(1), LARGE(2) / O077777, O177777 /
C     DATA LARGE(3), LARGE(4) / O177777, O177777 /
C     DATA RIGHT(1), RIGHT(2) / O022200, O000000 /
C     DATA RIGHT(3), RIGHT(4) / O000000, O000000 /
C     DATA DIVER(1), DIVER(2) / O022400, O000000 /
C     DATA DIVER(3), DIVER(4) / O000000, O000000 /
C     DATA LOG10(1), LOG10(2) / O037632, O020232 /
C     DATA LOG10(3), LOG10(4) / O102373, O147770 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE SUN
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE SUN
C     USING THE -r8 COMPILER OPTION
C
C     DATA DMACH(1) / Z'00010000000000000000000000000000' /
C     DATA DMACH(2) / Z'7FFEFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3F8E0000000000000000000000000000' /
C     DATA DMACH(4) / Z'3F8F0000000000000000000000000000' /
C     DATA DMACH(5) / Z'3FFD34413509F79FEF311F12B35816F9' /
C
C     MACHINE CONSTANTS FOR THE SUN 386i
C
C     DATA SMALL(1), SMALL(2) / Z'FFFFFFFD', Z'000FFFFF' /
C     DATA LARGE(1), LARGE(2) / Z'FFFFFFB0', Z'7FEFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'000000B0', Z'3CA00000' /
C     DATA DIVER(1), DIVER(2) / Z'FFFFFFCB', Z'3CAFFFFF'
C     DATA LOG10(1), LOG10(2) / Z'509F79E9', Z'3FD34413' /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
C
C     DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 /
C     DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 /
C     DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 /
C     DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 /
C     DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 /
C
C***FIRST EXECUTABLE STATEMENT  D1MACH
      IF (I .LT. 1 .OR. I .GT. 5) CALL XERMSG ('SLATEC', 'D1MACH',
     +   'I OUT OF BOUNDS', 1, 2)
C
      D1MACH = DMACH(I)
      RETURN
C
      END
      DOUBLE PRECISION FUNCTION DGAMLN(Z,IERR)
C***BEGIN PROLOGUE  DGAMLN
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  830501   (YYMMDD)
C***CATEGORY NO.  B5F
C***KEYWORDS  GAMMA FUNCTION,LOGARITHM OF GAMMA FUNCTION
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE LOGARITHM OF THE GAMMA FUNCTION
C***DESCRIPTION
C
C               **** A DOUBLE PRECISION ROUTINE ****
C         DGAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR
C         Z.GT.0.  THE ASYMPTOTIC EXPANSION IS USED TO GENERATE VALUES
C         GREATER THAN ZMIN WHICH ARE ADJUSTED BY THE RECURSION
C         G(Z+1)=Z*G(Z) FOR Z.LE.ZMIN.  THE FUNCTION WAS MADE AS
C         PORTABLE AS POSSIBLE BY COMPUTIMG ZMIN FROM THE NUMBER OF BASE
C         10 DIGITS IN A WORD, RLN=AMAX1(-ALOG10(R1MACH(4)),0.5E-18)
C         LIMITED TO 18 DIGITS OF (RELATIVE) ACCURACY.
C
C         SINCE INTEGER ARGUMENTS ARE COMMON, A TABLE LOOK UP ON 100
C         VALUES IS USED FOR SPEED OF EXECUTION.
C
C     DESCRIPTION OF ARGUMENTS
C
C         INPUT      Z IS D0UBLE PRECISION
C           Z      - ARGUMENT, Z.GT.0.0D0
C
C         OUTPUT      DGAMLN IS DOUBLE PRECISION
C           DGAMLN  - NATURAL LOG OF THE GAMMA FUNCTION AT Z.NE.0.0D0
C           IERR    - ERROR FLAG
C                     IERR=0, NORMAL RETURN, COMPUTATION COMPLETED
C                     IERR=1, Z.LE.0.0D0,    NO COMPUTATION
C
C
C***REFERENCES  COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C***ROUTINES CALLED  I1MACH,D1MACH
C***END PROLOGUE  DGAMLN
      DOUBLE PRECISION CF, CON, FLN, FZ, GLN, RLN, S, TLG, TRM, TST,
     * T1, WDTOL, Z, ZDMY, ZINC, ZM, ZMIN, ZP, ZSQ, D1MACH
      INTEGER I, IERR, I1M, K, MZ, NZ, I1MACH
      DIMENSION CF(22), GLN(100)
C           LNGAMMA(N), N=1,100
      DATA GLN(1), GLN(2), GLN(3), GLN(4), GLN(5), GLN(6), GLN(7),
     1     GLN(8), GLN(9), GLN(10), GLN(11), GLN(12), GLN(13), GLN(14),
     2     GLN(15), GLN(16), GLN(17), GLN(18), GLN(19), GLN(20),
     3     GLN(21), GLN(22)/
     4     0.00000000000000000D+00,     0.00000000000000000D+00,
     5     6.93147180559945309D-01,     1.79175946922805500D+00,
     6     3.17805383034794562D+00,     4.78749174278204599D+00,
     7     6.57925121201010100D+00,     8.52516136106541430D+00,
     8     1.06046029027452502D+01,     1.28018274800814696D+01,
     9     1.51044125730755153D+01,     1.75023078458738858D+01,
     A     1.99872144956618861D+01,     2.25521638531234229D+01,
     B     2.51912211827386815D+01,     2.78992713838408916D+01,
     C     3.06718601060806728D+01,     3.35050734501368889D+01,
     D     3.63954452080330536D+01,     3.93398841871994940D+01,
     E     4.23356164607534850D+01,     4.53801388984769080D+01/
      DATA GLN(23), GLN(24), GLN(25), GLN(26), GLN(27), GLN(28),
     1     GLN(29), GLN(30), GLN(31), GLN(32), GLN(33), GLN(34),
     2     GLN(35), GLN(36), GLN(37), GLN(38), GLN(39), GLN(40),
     3     GLN(41), GLN(42), GLN(43), GLN(44)/
     4     4.84711813518352239D+01,     5.16066755677643736D+01,
     5     5.47847293981123192D+01,     5.80036052229805199D+01,
     6     6.12617017610020020D+01,     6.45575386270063311D+01,
     7     6.78897431371815350D+01,     7.12570389671680090D+01,
     8     7.46582363488301644D+01,     7.80922235533153106D+01,
     9     8.15579594561150372D+01,     8.50544670175815174D+01,
     A     8.85808275421976788D+01,     9.21361756036870925D+01,
     B     9.57196945421432025D+01,     9.93306124547874269D+01,
     C     1.02968198614513813D+02,     1.06631760260643459D+02,
     D     1.10320639714757395D+02,     1.14034211781461703D+02,
     E     1.17771881399745072D+02,     1.21533081515438634D+02/
      DATA GLN(45), GLN(46), GLN(47), GLN(48), GLN(49), GLN(50),
     1     GLN(51), GLN(52), GLN(53), GLN(54), GLN(55), GLN(56),
     2     GLN(57), GLN(58), GLN(59), GLN(60), GLN(61), GLN(62),
     3     GLN(63), GLN(64), GLN(65), GLN(66)/
     4     1.25317271149356895D+02,     1.29123933639127215D+02,
     5     1.32952575035616310D+02,     1.36802722637326368D+02,
     6     1.40673923648234259D+02,     1.44565743946344886D+02,
     7     1.48477766951773032D+02,     1.52409592584497358D+02,
     8     1.56360836303078785D+02,     1.60331128216630907D+02,
     9     1.64320112263195181D+02,     1.68327445448427652D+02,
     A     1.72352797139162802D+02,     1.76395848406997352D+02,
     B     1.80456291417543771D+02,     1.84533828861449491D+02,
     C     1.88628173423671591D+02,     1.92739047287844902D+02,
     D     1.96866181672889994D+02,     2.01009316399281527D+02,
     E     2.05168199482641199D+02,     2.09342586752536836D+02/
      DATA GLN(67), GLN(68), GLN(69), GLN(70), GLN(71), GLN(72),
     1     GLN(73), GLN(74), GLN(75), GLN(76), GLN(77), GLN(78),
     2     GLN(79), GLN(80), GLN(81), GLN(82), GLN(83), GLN(84),
     3     GLN(85), GLN(86), GLN(87), GLN(88)/
     4     2.13532241494563261D+02,     2.17736934113954227D+02,
     5     2.21956441819130334D+02,     2.26190548323727593D+02,
     6     2.30439043565776952D+02,     2.34701723442818268D+02,
     7     2.38978389561834323D+02,     2.43268849002982714D+02,
     8     2.47572914096186884D+02,     2.51890402209723194D+02,
     9     2.56221135550009525D+02,     2.60564940971863209D+02,
     A     2.64921649798552801D+02,     2.69291097651019823D+02,
     B     2.73673124285693704D+02,     2.78067573440366143D+02,
     C     2.82474292687630396D+02,     2.86893133295426994D+02,
     D     2.91323950094270308D+02,     2.95766601350760624D+02,
     E     3.00220948647014132D+02,     3.04686856765668715D+02/
      DATA GLN(89), GLN(90), GLN(91), GLN(92), GLN(93), GLN(94),
     1     GLN(95), GLN(96), GLN(97), GLN(98), GLN(99), GLN(100)/
     2     3.09164193580146922D+02,     3.13652829949879062D+02,
     3     3.18152639620209327D+02,     3.22663499126726177D+02,
     4     3.27185287703775217D+02,     3.31717887196928473D+02,
     5     3.36261181979198477D+02,     3.40815058870799018D+02,
     6     3.45379407062266854D+02,     3.49954118040770237D+02,
     7     3.54539085519440809D+02,     3.59134205369575399D+02/
C             COEFFICIENTS OF ASYMPTOTIC EXPANSION
      DATA CF(1), CF(2), CF(3), CF(4), CF(5), CF(6), CF(7), CF(8),
     1     CF(9), CF(10), CF(11), CF(12), CF(13), CF(14), CF(15),
     2     CF(16), CF(17), CF(18), CF(19), CF(20), CF(21), CF(22)/
     3     8.33333333333333333D-02,    -2.77777777777777778D-03,
     4     7.93650793650793651D-04,    -5.95238095238095238D-04,
     5     8.41750841750841751D-04,    -1.91752691752691753D-03,
     6     6.41025641025641026D-03,    -2.95506535947712418D-02,
     7     1.79644372368830573D-01,    -1.39243221690590112D+00,
     8     1.34028640441683920D+01,    -1.56848284626002017D+02,
     9     2.19310333333333333D+03,    -3.61087712537249894D+04,
     A     6.91472268851313067D+05,    -1.52382215394074162D+07,
     B     3.82900751391414141D+08,    -1.08822660357843911D+10,
     C     3.47320283765002252D+11,    -1.23696021422692745D+13,
     D     4.88788064793079335D+14,    -2.13203339609193739D+16/
C
C             LN(2*PI)
      DATA CON                    /     1.83787706640934548D+00/
C
C***FIRST EXECUTABLE STATEMENT  DGAMLN
      IERR=0
      IF (Z.LE.0.0D0) GO TO 70
      IF (Z.GT.101.0D0) GO TO 10
      NZ = INT(SNGL(Z))
      FZ = Z - FLOAT(NZ)
      IF (FZ.GT.0.0D0) GO TO 10
      IF (NZ.GT.100) GO TO 10
      DGAMLN = GLN(NZ)
      RETURN
   10 CONTINUE
      WDTOL = D1MACH(4)
      WDTOL = DMAX1(WDTOL,0.5D-18)
      I1M = I1MACH(14)
      RLN = D1MACH(5)*FLOAT(I1M)
      FLN = DMIN1(RLN,20.0D0)
      FLN = DMAX1(FLN,3.0D0)
      FLN = FLN - 3.0D0
      ZM = 1.8000D0 + 0.3875D0*FLN
      MZ = INT(SNGL(ZM)) + 1
      ZMIN = FLOAT(MZ)
      ZDMY = Z
      ZINC = 0.0D0
      IF (Z.GE.ZMIN) GO TO 20
      ZINC = ZMIN - FLOAT(NZ)
      ZDMY = Z + ZINC
   20 CONTINUE
      ZP = 1.0D0/ZDMY
      T1 = CF(1)*ZP
      S = T1
      IF (ZP.LT.WDTOL) GO TO 40
      ZSQ = ZP*ZP
      TST = T1*WDTOL
      DO 30 K=2,22
        ZP = ZP*ZSQ
        TRM = CF(K)*ZP
        IF (DABS(TRM).LT.TST) GO TO 40
        S = S + TRM
   30 CONTINUE
   40 CONTINUE
      IF (ZINC.NE.0.0D0) GO TO 50
      TLG = DLOG(Z)
      DGAMLN = Z*(TLG-1.0D0) + 0.5D0*(CON-TLG) + S
      RETURN
   50 CONTINUE
      ZP = 1.0D0
      NZ = INT(SNGL(ZINC))
      DO 60 I=1,NZ
        ZP = ZP*(Z+FLOAT(I-1))
   60 CONTINUE
      TLG = DLOG(ZDMY)
      DGAMLN = ZDMY*(TLG-1.0D0) - DLOG(ZP) + 0.5D0*(CON-TLG) + S
      RETURN
C
C
   70 CONTINUE
      IERR=1
      RETURN
      END
      SUBROUTINE DQG16(XL,XU,FCT,Y)
      IMPLICIT NONE
C
C SUBROUTINE DQG16 COMPUTES INTEGRAL(FCT(X),SUMMED OVER X FROM XL TO XU)
C USING 16-POINT GAUSSIAN QUADRATURE FORMULA, WHICH INTEGRATES
C POLYNOMIALS OF UP TO DEGREE 31 EXACTLY
C THIS SUBROUTINE IS COPIED FROM IBM SSP
C USAGE: CALL DQG16(XL,XU,FCT,Y)
C        PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT IN CALLING PROGRAM
C        XL=DOUBLE PRECISION LOWER BOUND OF INTERVAL
C        XU=DOUBLE PRECISION UPPER BOUND OF INTERVAL
C        FCT=NAME OF EXTERNAL DOUBLE PRECISION FUNCTION
C        Y=RESULTING DOUBLE PRECISION INTEGRAL VALUE
C
      DOUBLE PRECISION XL,XU,Y,A,B,C,FCT
      EXTERNAL FCT
      A=.5D0*(XU+XL)
      B=XU-XL
      C=.49470046749582497D0*B
      Y=.13576229705877047D-1*(FCT(A+C)+FCT(A-C))
      C=.47228751153661629D0*B
      Y=Y+.31126761969323946D-1*(FCT(A+C)+FCT(A-C))
      C=.43281560119391587D0*B
      Y=Y+.47579255841246392D-1*(FCT(A+C)+FCT(A-C))
      C=.37770220417750152D0*B
      Y=Y+.62314485627766936D-1*(FCT(A+C)+FCT(A-C))
      C=.30893812220132187D0*B
      Y=Y+.7479799440828837D-1*(FCT(A+C)+FCT(A-C))
      C=.22900838882861369D0*B
      Y=Y+.8457825969750127D-1*(FCT(A+C)+FCT(A-C))
      C=.14080177538962946D0*B
      Y=Y+.9130170752246179D-1*(FCT(A+C)+FCT(A-C))
      C=.47506254918818720D-1*B
      Y=B*(Y+.9472530522753425D-1*(FCT(A+C)+FCT(A-C)))
      RETURN
      END  
      SUBROUTINE EFFMED(ITHEORY,FINC,EPSMAT,EPSINC,EPSEFF)
      IMPLICIT NONE

c arguments

      INTEGER ITHEORY
      DOUBLE PRECISION FINC
      DOUBLE COMPLEX EPSEFF,EPSINC,EPSMAT

c local variables

      DOUBLE COMPLEX B
c-----------------------------------------------------------------------
c subroutine EFFMED
c given:
c
c     ITHEORY = 0 for simple volume-weighted average
c               1 for Bruggemann theory
c               2 for Maxwell-Garnett theory
c     FINC   = volume filling factor of inclusions
c     EPSMAT = complex dielectric function for matrix material
c     EPSINC = complex dielectric function for inclusions
c
c returns
c
c     EPSEFF = complex dielectric function for effective medium
c
C B.T.Draine, Princeton Univ. Observatory, 94.04.29
C history
C 02.12.22 (BTD) converted to subroutine
c end history
c-----------------------------------------------------------------------
c sanity check:

      IF(FINC.LT.0.D0.OR.FINC.GT.1.D0)THEN
         WRITE(0,*)'Fatal error in effmed: FINC=',FINC
         STOP
      ENDIF
      IF(ITHEORY.EQ.0)THEN

c simple linear average of dielectric functions:

         EPSEFF=(1.-FINC)*EPSMAT+FINC*EPSINC

      ELSEIF(ITHEORY.EQ.1)THEN

c Bruggemann theory:

         B=(1.D0-1.5D0*FINC)*EPSMAT+(1.5D0*FINC-.5D0)*EPSINC
         EPSEFF=.5D0*(B+SQRT(B*B+2.D0*EPSINC*EPSMAT))

      ELSEIF(ITHEORY.EQ.2)THEN

c Maxwell Garnett theory:

         EPSEFF=EPSMAT*(1.D0+3.D0*FINC*(EPSINC-EPSMAT)/
     &          (EPSINC+2.*EPSMAT-FINC*(EPSINC-EPSMAT)))
      ELSE

         WRITE(0,*)'Fatal error in effmed: option',
     &             ' ITHEORY=',ITHEORY,' not recognized'

         STOP
      ENDIF
      RETURN
      END
   
      SUBROUTINE EMDIPL(X,CXINDEX,CXMU,QABS,QSCA,QEXT)
      IMPLICIT NONE

c Arguments:

      REAL QABS,QEXT,QSCA,X
      COMPLEX CXINDEX,CXMU

c Local variables

      REAL MUFAC
      COMPLEX CXTERM
      DOUBLE PRECISION DX
      DOUBLE COMPLEX CXI,CXINDEX2,Y

      INTRINSIC DCMPLX,DIMAG

c***********************************************************************
c Given:
c     X = 2*pi*a/lambda
c     CXINDEX = complex refractive index
c

c Returns:
c     QABS = C_abs/(pi*a**2)
c     QSCA = C_sca/(pi*a**2)
c     QEXT = QABS+QSCA
c
c Using electric and magnetic dipole approximation for a sphere
c consisting of nonmagnetic material.
c
c B.T. Draine, Princeton University Observatory
c History:
c Created sometime in 1980s
c 97.02.18 (BTD) Added SAVE statement to make f77 compliant and
c                g77 compatible
c 99.11.01 (BTD) groomed to please FTNCHEK
c 11.03.03 (BTD) v2
c                * add CXMU to argument list
c                * add effect of CXMU on eddy current dissipation
c 11.03.05 (BTD) * corrected bugs                 
c 11.12.12 (BTD) * edited comments
c end history
c***********************************************************************
      CXI=(0.D0,1.D0)
      DX=DBLE(X)
      Y=DCMPLX(CXINDEX*X)

c Series expansion in X for small X

      CXINDEX2=DCMPLX(CXINDEX*CXINDEX)

c Term proportional to X is electric dipole absorption:

      QABS=REAL(4.D0*DX*DIMAG((CXINDEX2-1.D0)/(CXINDEX2+2.D0)))

c QABSM is the contribution to QABS from magnetic dipole absorption
c*****************************************************************
c Note: van de Hulst (1957), section 14.21 has given X**3 term in
c expansion of Mie cross section in powers of X.  We do NOT use this
c as it is not a particularly good estimate if one is dealing with
c highly conductive particles with Y > 1 even when X << 1 ???
c      QABSM=-(4.E0/15.E0)*X**3*
c     &       DIMAG(((CXINDEX2-1.E0)/(CXINDEX2+2.E0))**2*
c     &       (38.E0+CXINDEX2*(27.E0+CXINDEX2))/(2.E0*CXINDEX2+3.E0))
c******************************************************************
c Our estimate for QABSM is not exact to order X**3 (it does not
c include correction to electric dipole absorption) but does
c accurately represent magnetic dipole absorption even when |m|x > 1.
c Result for magnetic polarizability of conducting sphere is taken
c from Landau & Lifshitz, Electrodynamics of Continuous Media,
c sections 45, 72, and 73.

c Formula in Landau & Liftshitz for magnetic dipole dissipation
c is for nonmagnetic materials.  Dissipation is due to joule heating
c by eddy currents induced by time-varying B field.
c Applied EM wave drives eddy currents according to curl E = -(1/c)dB/dt
c In nonmagnetic medium, B = H
c In magnetic particle, B is enhanced by magnetic susceptibility
c In magnetostatics limit, H within the sphere is related to the applied
c H_o by
c        1        H_o
c H_in = - * ------------  where L=1/3 for a sphere
c        L   mu - 1 + 1/L
c
c         3 H_o
c H_in = ------
c        mu + 2
c
c        3*mu*H_0
c B_in = --------
c         mu + 2
c
c The induced E [from curl E = -(1/c)dB/dt] should be
c proportional to omega*B_in
c
c The current J should be proportional to E
c Hence the joule heating propto J^2
c should be proportional to |B_in|^2

c We therefore take LL formula for magnetic dipole absorption and
c multiply it by |3*mu/(2+mu)|^2

c However, the magnetic permeability also contributes to the
c "backreaction" term that would limit the dissipation if
c sigma -> infinity or mu -> infinity

      CXTERM=3.*CXMU/(CXMU+2.)

      MUFAC=REAL(CXTERM*CONJG(CXTERM))

c*** diagnostic
c      write(0,*)'emdipl_v2 ckpt 10'
c      write(0,*)' cxmu=',cxmu
c      write(0,*)' mufac=',mufac
c**
      IF(ABS(Y).LT.1.D-1)THEN

c Expansion for small Y:

         QABS=QABS+MUFAC*REAL(2.D0*DX*DIMAG(Y*Y)/15.D0)
      ELSE
         IF(DIMAG(Y).LT.5.D0)THEN

c General formula:

            QABS=QABS+
     &           MUFAC*REAL(6.D0*DX*DIMAG((1.D0/Y-COS(Y)/SIN(Y))/Y))

         ELSE

c Approximation for Im(Y)>5: cot(Y)=-I to high accuracy

            QABS=QABS+MUFAC*REAL(6.D0*DX*DIMAG((1.D0/Y+CXI)/Y))

         ENDIF
      ENDIF

c Magnetic dipole scattering:
c Landau & Lifshitz formula for the magnetic moment is for the magnetic
c dipole arising from induced currents.
c In addition there will be magnetic moment from aligned spins.
c The ???

      QSCA=REAL((8.D0/3.D0)*ABS(((CXINDEX2-1.D0)/
     &          (CXINDEX2+2.D0))**2)*DX**4)
      QEXT=QSCA+QABS
      RETURN
      END
      SUBROUTINE ErrMsg( MESSAG, FATAL )

c        Print out a warning or error message;  abort if error
c        after making symbolic dump (machine-specific)

c     .. Scalar Arguments ..

      CHARACTER MESSAG*(*)
      LOGICAL   FATAL
c     ..
c     .. Local Scalars ..

      LOGICAL   MSGLIM
      INTEGER   MAXMSG, NUMMSG
c     ..
c     .. External Subroutines ..

ccccc EXTERNAL  SYMDUMP
c     ..
      SAVE      MAXMSG, NUMMSG, MSGLIM
c      DATA      NUMMSG / 0 /,  MAXMSG / 100 /,  MSGLIM /.FALSE./

      DATA NUMMSG/0/,MAXMSG/2/,MSGLIM/.FALSE./


      IF( FATAL ) THEN

         WRITE(0, '(//,2A,//)' ) ' ****** ERROR *****  ', MESSAG

c                                 ** Example symbolic dump call for Cray
ccccc    CALL SYMDUMP( '-B -c3' )

         STOP

      END IF


      NUMMSG = NUMMSG + 1

      IF( MSGLIM ) RETURN

      IF( NUMMSG.LE.MAXMSG ) THEN

c         WRITE(0, '(/,2A,/)' ) ' ****** WARNING *****  ', MESSAG
         WRITE(0,'(2A)')' ****** WARNING *****  ',MESSAG

      ELSE

         WRITE(0, 9000 )
         MSGLIM = .True.
      END IF


      RETURN

c 9000 FORMAT( / , / , ' ****** TOO MANY WARNING MESSAGES --  ',
c     &      'They will no longer be printed *******', / , / )
 9000 FORMAT(' ****** TOO MANY WARNING MESSAGES --  ',
     &      'They will no longer be printed *******')
      END

      LOGICAL FUNCTION WrtBad( VarNam )

c          Write names of erroneous variables and return 'TRUE'

c      INPUT :   VarNam = Name of erroneous variable to be written
c                         ( CHARACTER, any length )

c     .. Scalar Arguments ..

      CHARACTER VarNam*(*)
c     ..
c     .. Local Scalars ..

      INTEGER   MAXMSG, NUMMSG
c     ..
c     .. External Subroutines ..

      EXTERNAL  ErrMsg
c     ..
      SAVE      NUMMSG, MAXMSG
      DATA      NUMMSG / 0 /, MAXMSG / 50 /


      WrtBad = .TRUE.
      NUMMSG = NUMMSG + 1
      WRITE(0, '(3A)' ) ' ****  Input variable  ', VarNam,
     &                   '  in error  ****'

      IF( NUMMSG.EQ.MAXMSG )
     &    CALL ErrMsg( 'Too many input errors.  Aborting...',.TRUE.)

      RETURN
      END

      LOGICAL FUNCTION WrtDim( DimNam, Minval )

c          Write name of too-small symbolic dimension and
c          the value it should be increased to;  return 'TRUE'

c      INPUT :  DimNam = Name of symbolic dimension which is too small
c                        ( CHARACTER, any length )
c               Minval = Value to which that dimension should be
c                        increased (at least)

c     .. Scalar Arguments ..

      CHARACTER DimNam*(*)
      INTEGER   Minval
c     ..

      WRITE(0, '(3A,I7)' ) ' ****  Symbolic dimension  ', DimNam,
     &                      '  should be increased to at least ', Minval
      WrtDim = .TRUE.
      RETURN
      END

      LOGICAL FUNCTION TstBad( VarNam, RelErr )

c       Write name (VarNam) of variable failing self-test and its
c       percent error from the correct value;  return  'FALSE'.

c     .. Scalar Arguments ..

      CHARACTER VarNam*(*)
c      REAL      RelErr
      DOUBLE PRECISION RelErr
c     ..

      TstBad = .FALSE.
      WRITE(0, '(/,3A,1P,E11.2,A)' ) ' *** Output variable ', VarNam,
     &   ' differed by ', 100.*RelErr,
     &   ' per cent from correct value.  Self-test failed.'
      RETURN
      END
*DECK FDUMP
      SUBROUTINE FDUMP
C***BEGIN PROLOGUE  FDUMP
C***PURPOSE  Symbolic dump (should be locally written).
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3
C***TYPE      ALL (FDUMP-A)
C***KEYWORDS  ERROR, XERMSG
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C        ***Note*** Machine Dependent Routine
C        FDUMP is intended to be replaced by a locally written
C        version which produces a symbolic dump.  Failing this,
C        it should be replaced by a version which prints the
C        subprogram nesting list.  Note that this dump must be
C        printed on each of up to five files, as indicated by the
C        XGETUA routine.  See XSETUA and XGETUA for details.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  FDUMP
C***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
      SUBROUTINE FINDPATH(PATH)
      IMPLICIT NONE
      CHARACTER*15 PATH
      PATH='/u/draine/work/'
      RETURN
      END
      SUBROUTINE FREE(EEV,IDIEL,TGR,AGR,DEPS1,DEPS2)
      IMPLICIT NONE
c---------------------- subroutine free v4 ----------------------------
C Arguments:
      INTEGER IDIEL
      REAL AGR,DEPS1,DEPS2,EEV,TGR
C Local variables:
      LOGICAL INIT
      INTEGER NRESMX
      PARAMETER(NRESMX=13)
      INTEGER JR,NRES
      REAL EPL,FOURPI,GAMMA,OMEGA,SIGDC,TAU,TAUB,TERM,VF,X,Y
      REAL 
     &   ERES(1:NRESMX),
     &   GAMM(1:NRESMX),
     &   SRES(1:NRESMX)
      COMPLEX CXI,DCXEPS
C***********************************************************************
C Subroutine FREE computes contributions of free electrons to
C dielectric function of graphite or other materials with free-electron
c contribution.
c AND
c adds contributions of specified resonances
c
C Returns DEPS1=DEPS2=0 when called for silicate, diamond, SiC.
C Input:
C     EEV  = photon energy/eV
C     IDIEL=1 for graphite, E parallel to c axis
C             DL84 conductivity sigma_dc=29.1 mho/cm
C           2 for graphite, E perpendicular to c axis, DL84 version
C           9 for new graphite, E parallel to c axis, with
C                 increased dc conductivity sigma_dc=62.5 mho/cm
C          12 for D03 graphite E parallel to c axis, 
C                 increased dc conductivity sigma_dc=62.5 mho/cm
C          13 for D03 graphite E perpendicular to c axis, 
C                 DL84 conductivity
C          15 for D03 graphite, E parallel to c, with porosity(a)
C                 sigma_dc=62.5 mho/cm
C                 porosity effects added later in subroutine INDEX
C          16 for D03 graphite, E perp. to c, with porosity(a)
C                 porosity effects added later in subroutine INDEX
C          17 for Cpaexp1: E parallel to c. Djuric & Li 99, optical data
C                 free electron contribution with sigma_dc(298K)=62.5 mho/cm
C          18 for Cpaexp2: E parallel to c. Djuric & Li 99, EELS data
C                 free electron contribution with sigma_dc(298K)=62.5 mho/cm
C          19 for Cpeexp: E perp to c. Djuric & Li 99
C                 free electron contribution with sigma_dc(298K)=7830 mho/cm
C
C          29 for amcACH2: modified ACH2 dielectric function from
C                          Zubko et al
C                          need to add free electron contribution
C          30 for graphite, E||c, with new two-component representation
C                           for the free-electron contribution
C          31 for Fe metal
c          32 for Fe3O4 magnetite
C
C     TGR  = grain temperature (K)
C     AGR  = grain radius (micron)
C Returns:
C     DEPS1=contribution to Re(epsilon) due to free electrons
C     DEPS2=contribution to Im(epsilon) due to free electrons
C
C Free electrons are assumed to be described by Drude model.
C B.T.Draine, Princeton Univ. Obs., 83.09
C History:
C 83.11.08 (BTD) modified.
C 86.08.11 (BTD) change tau(parallel) from 1.e-14 to 1.4e-14 sec
C 92.01.06 (BTD) explicit declaration of all variables.
C 00.07.18 (BTD) modified to support new silicate experiments (IDIEL=8)
C 00.09.13 (BTD) modified to support new graphite experiments (IDIEL=9)
C 00.09.14 (BTD) correction to support IDIEL=9
C 00.10.27 (BTD) cosmetic change 
C 02.12.22 (BTD) modified to support IDIEL=18,19,21,22
C 03.01.03 (BTD) modified to change IDIEL identifications, now support
C                IDIEL=17-19, 22-24, 31-33
C 07.11.16 (BTD) modified to recognize IDIEL=43 (amcACH2)
C                and provide free-electron contribution
c 07.12.09 (BTD) had been using TAU=3.e-14 for IDIEL=9,18,23
c                which had been intended as an *experiment*
c                but was inadvertently left in code
c                changed back to 1.4e-14 on 07.12.09
c 07.12.22 (BTD) Upon consideration,decided that we actually
c                want to keep TAU=3.e-14 for IDIEL=9
c 11.02.13 (BTD) v3
c                * Add Fe metal: 3 component model for free electrons
c 11.02.17 (BTD) * revise parameters for Fe
c 12.05.21 (BTD) v4
c                * add Fe3O4 magnetite free-electron term
C 12.05.22 (BTD) * modify Fe3O4 magnetite free-electron term
c                * add resonances for Fe3O4
c                * add initialization section for CXI,FOURPI
c 12.05.24 (BTD) * modify sigma_dc for magnetite
C end history
C***********************************************************************
      SAVE CXI,FOURPI,INIT
      DATA INIT/.TRUE./
      IF(INIT)THEN
         CXI=(0.E0,1.E0)
         FOURPI=16.D0*ATAN(1.D0)
         INIT=.FALSE.
      ENDIF
      DEPS1=0.E0
      DEPS2=0.E0

C nonzero contribution to DEPS1,DEPS2 only for following cases:
C   IDIEL = 1 DL84 graphite, E parallel to c axis
C           2 DL84 graphite, E perp to c axis
C           9 experimental modification to graphite, E parallel to c
C          18 D03 graphite, E parallel to c
C          19 D03 graphite, E perp to c
C          23 D03 graphite, E parallel to c
C          24 D03 graphite, E perp to c
c          30
c          31 Fe metal
c          32 magnetite Fe3O4

C Do we need to add a free-electron contribution?

      IF(IDIEL.NE.1.AND.
     &   IDIEL.NE.2.AND.
     &   IDIEL.NE.9.AND.
     &   IDIEL.NE.12.AND.
     &   IDIEL.NE.13.AND.
     &   IDIEL.NE.15.AND.
     &   IDIEL.NE.16.AND.
     &   IDIEL.NE.17.AND.
     &   IDIEL.NE.18.AND.
     &   IDIEL.NE.19.AND.
     &   IDIEL.NE.29.AND.
     &   IDIEL.NE.30.AND.
     &   IDIEL.NE.31.AND.
     &   IDIEL.NE.32)RETURN

c----- various versions of graphite, E || c -----

      IF(IDIEL.EQ.1)THEN

c DL84 graphite, E || c
c omega_p = .101 eV = 1.5345e14 s-1 [independent of TGR]
c tau_bulk = 1.4e-14 s-1
c bulk conductivity sigma_0 = omega_p^2*tau_bulk/4*pi=2.6232e13 s-1
c      = (2.6232e13/9e11) mho cm-1
c      = 29.1 mho cm-1
c For motion along c axis, VF taken to be 3.7e6 cm/s=3.7e10 micron/s
c at low temperatures, based on assumed Fermi energy=0.022eV
c (Williamson, Foner, and Dresselhaus 1965, Phys.Rev.A140,1429)
c and electron effective mass=5.7m_e (Soule, McClure, and Smith 1964,
c Phys.Rev.A134,453)

         VF=3.7E10*SQRT(1.E0+TGR/2.55E2)

         EPL=0.101
         TAUB=1.4E-14
         TAU=1.E0/(1.E0/TAUB+VF/AGR)
         Y=EEV/EPL
         GAMMA=1.E0/(1.518E15*EPL*TAU)
         DEPS1=-1.E0/(Y*Y+GAMMA*GAMMA)
         DEPS2=GAMMA/(Y*(Y*Y+GAMMA*GAMMA))
         RETURN

      ELSEIF(IDIEL.EQ.2.OR.
     &       IDIEL.EQ.13.OR.
     &       IDIEL.EQ.16.OR.
     &       IDIEL.EQ.19)THEN

c graphite, E perp to c, free electron contribution following DL84.
c EPL is fit by a quadratic following DL84:
c at 298K, EPL=0.4407 eV, omega_p=6.692e14 s-1
         
         EPL=.285E0*SQRT(1.E0-6.24E-3*TGR+3.66E-5*TGR*TGR)

c 1/tau_bulk is fit by a quadratic, following DL84
c at 298K, tau_bulk=1.977e-13 s

         TAUB=4.2E-11/(1.E0+.322E0*TGR+1.30E-3*TGR*TGR)

c From DL84:
c At 298K, sigma_dc = omega_p^2*tau_bulk/4*pi=7.045e15 s-1
c                   = 7.828e3 mho/cm  
c [experimental value = 2.5e4 mho/cm (Soule 1958)]
c From DL84:
c For motion in basal plane, VF taken to be 4.5e7 cm/s=4.5e11 micron/s
c at low temperature, based on assumed Fermi energy 0.022 eV
c (Williamson, Foner, and Dresselhaus 1965, Phys.Rev.A140,1429)
c and electron effective mass=0.039m_e (Soule, McClure, and Smith 1964,
c Phys.Rev.A134,453).
c Include correction for thermal excitation: 255 degK=.022eV/k

         VF=4.5E11*SQRT(1.E0+TGR/2.55E2)

         TAU=1.E0/(1.E0/TAUB+VF/AGR)
         Y=EEV/EPL
         GAMMA=1.E0/(1.518E15*EPL*TAU)
         DEPS1=DEPS1-1.E0/(Y*Y+GAMMA*GAMMA)
         DEPS2=DEPS2+GAMMA/(Y*(Y*Y+GAMMA*GAMMA))
         RETURN

      ELSEIF(IDIEL.EQ.9.OR.
     &       IDIEL.EQ.12.OR.
     &       IDIEL.EQ.15.OR.
     &       IDIEL.EQ.17.OR.
     &       IDIEL.EQ.18)THEN

c graphite, E parallel to c, enhanced dc conductivity
c
c for 9,12,15,17,18: experiment with new value of tau_bulk:
c                    tau_bulk increased to 3.e-14 s
c in order to increase d.c. conductivity by factor 3/1.4
c              sigma_dc = 5.6211e13 s-1 = 62.5 mho cm-1
c This d.c. conductivity is still within the range of experimental
c values cited by DL84.
                        
c Objective is to **decrease** lambda > 100um emissivity
c in order to improve agreement with observed lambda > 100um emission
c spectrum of high latitude MW dust observed by Finkbeiner, Davis,
c & Schlegel 1999 [see Fig. 14 of Draine & Li 2007, which was computed
c for tau_bulk=3.e-14]

c free electron velocity, in micron/s

         VF=3.7E10*SQRT(1.E0+TGR/2.55E2)

         EPL=0.101
         TAUB=3.0E-14
         TAU=1.E0/(1.E0/TAUB+VF/AGR)
         Y=EEV/EPL
         GAMMA=1.E0/(1.518E15*EPL*TAU)
         DEPS1=-1.E0/(Y*Y+GAMMA*GAMMA)
         DEPS2=GAMMA/(Y*(Y*Y+GAMMA*GAMMA))
         RETURN

      ELSEIF(IDIEL.EQ.29)THEN

c amorphous carbon amcACH2, based on ACH2 from Zubko et al 1996:
c we adopt sigma_0(mks)=16.3 ohm^{-1} cm^{-1}
c          sigma_0(cgs)=1.465e13 s-1
c we adopt tau=1e-14 s
c     ->   omega_p = (4*pi*sigma_0/tau)^{1/2} = 1.357e14 s-1
c          hbar*omega_p = 0.08935 eV

c take a wild guess at v_F=1.e7 cm/s=1.e11 um/s

         VF=1.E11
         TAUB=1.E-14
         EPL=0.08935

         TAU=1.E0/(1.E0/TAUB+VF/AGR)
         Y=EEV/EPL
         GAMMA=1.E0/(1.518E15*EPL*TAU)
         DEPS1=DEPS1-1.E0/(Y*Y+GAMMA*GAMMA)
         DEPS2=DEPS2+GAMMA/(Y*(Y*Y+GAMMA*GAMMA))
         RETURN

      ELSEIF(IDIEL.EQ.30)THEN

c modification to free electron component to reproduce
c d.c. conductivity of IDIEL=9
c             sigma_dc = 5.6211e13 s-1 = 62.5 mho cm-1
c but **suppress** 33um absorption peak present in IDIEL=9 and 18
c [which uses omega_p = 1.53e14 s-1 and tau_bulk = 3.e-14 s]
c this is accomplished by using two components for the free-electron
c contribution.
c each component is assumed to contribute 50% of the d.c. conductivity
c comp 1: omega_p = 1.e14 s-1  or  E_p = 0.0660 eV
c        tau_bulk = 0.5*(1.53/1)^2*3e-14
c                 = 3.5114e-14 s
c comp 2: omega_p = 2.e14 s-1  or  E_p = 0.132  eV
c        tau_bulk = 0.5*(1.53/2)^2*3e-14
c                 = 0.8778e-14 s

c free electron velocity, in micron/s

         VF=3.7E10*SQRT(1.E0+TGR/2.55E2)

         EPL=.0660
         TAUB=3.51E-14
         TAU=1.E0/(1.E0/TAUB+VF/AGR)
         Y=EEV/EPL
         GAMMA=1.E0/(1.518E15*EPL*TAU)
         DEPS1=-1.E0/(Y*Y+GAMMA*GAMMA)
         DEPS2=GAMMA/(Y*(Y*Y+GAMMA*GAMMA))
        
         EPL=0.132
         TAUB=8.778E-15
         TAU=1.E0/(1.E0/TAUB+VF/AGR)
         Y=EEV/EPL
         GAMMA=1.E0/(1.518E15*EPL*TAU)
         DEPS1=DEPS1-1.E0/(Y*Y+GAMMA*GAMMA)
         DEPS2=DEPS2+GAMMA/(Y*(Y*Y+GAMMA*GAMMA))

         RETURN

      ELSEIF(IDIEL.EQ.31)THEN

c Fe metal
c three-component model for free electron contribution in Fe
c guess v_F=1.98e8 cm/s = 1.98e12 um/s [Ashcroft & Mermin 1976]
c free electron model for bulk Fe (see /u/draine/work/opt/fe.sm):
c deps2 = 
c      hbar/tau_1    hbar*omega_p   tau_1 (s)
c  1    6.5 eV       14.42 eV       1.013e-16
c  2    0.0165 eV     3.146 eV      3.991e-14
c  3    0.010 eV      0.707 eV      6.586e-14
         VF=1.98E12

         TAUB=1.013E-16
         EPL=14.42
         TAU=1.E0/(1.E0/TAUB+VF/AGR)
         Y=EEV/EPL
         GAMMA=1.E0/(1.518E15*EPL*TAU)
         DEPS1=-1.E0/(Y*Y+GAMMA*GAMMA)
         DEPS2=GAMMA/(Y*(Y*Y+GAMMA*GAMMA))
c*** diagnostic
c         write(0,*)'free_v3 ckp 31.1: y=',y,' gamma=',gamma
c         write(0,*)'   deps1,deps2=',deps1,deps2
c***

         TAUB=3.991E-14
         EPL=3.146
         TAU=1.E0/(1.E0/TAUB+VF/AGR)
         Y=EEV/EPL
         GAMMA=1.E0/(1.518E15*EPL*TAU)
         DEPS1=DEPS1-1.E0/(Y*Y+GAMMA*GAMMA)
         DEPS2=DEPS2+GAMMA/(Y*(Y*Y+GAMMA*GAMMA))
c*** diagnostic
c         write(0,*)'free_v3 ckp 31.1: y=',y,' gamma=',gamma
c         write(0,*)'   deps1,deps2=',deps1,deps2
c***

         TAUB=6.59-14
         EPL=0.707
         TAU=1.E0/(1.E0/TAUB+VF/AGR)
         Y=EEV/EPL
         GAMMA=1.E0/(1.518E15*EPL*TAU)
         DEPS1=DEPS1-1.E0/(Y*Y+GAMMA*GAMMA)
         DEPS2=DEPS2+GAMMA/(Y*(Y*Y+GAMMA*GAMMA))
c*** diagnostic
c         write(0,*)'free_v3 ckp 31.1: y=',y,' gamma=',gamma
c         write(0,*)'   deps1,deps2=',deps1,deps2
c***

      ELSEIF(IDIEL.EQ.32)THEN

c Fe3O4 magnetite
c T_V = 119K = Verwey temperature
c Take the free-electron contribution to be
c sigma(dc)=3.e14*exp(-260/T) s-1 for T_V < T < 300K
c sigma(dc)=[3e16*exp(-1330/T)+3e12*exp(-740/T)] s-1 for T < T_V
c (Miles etal 1957; Tsuda et al 1991)
c
c adopt damping time tau_bulk=3.1e-13 sec  (hbar/tau = 0.00212 eV)
c for bulk material at T < 130K
c v_F = 1e8 cm/s = 1e12 um/s

         IF(TGR.LE.119.)THEN
            SIGDC=3.E12*EXP(-740./TGR)
            IF(TGR.GT.30.)SIGDC=SIGDC+3.E16*EXP(-1330/TGR)
         ELSEIF(TGR.GT.119..AND.TGR.LT.260.)THEN
            SIGDC=2.3E14*(TGR-100.)/160.
         ELSE
            SIGDC=2.3E14
         ENDIF
         OMEGA=1.518E15*EEV
         TAUB=3.1E-13
         TAU=1./(1./TAUB+1.E12/AGR)
         DCXEPS=FOURPI*CXI*(SIGDC*TAU/TAUB)/(OMEGA*(1.-CXI*OMEGA*TAU))

c*** diagnostic
c         write(0,fmt='(a,1pe10.3,a,1p2e10.3)')
c     &   'free_v4 ckpt 32.1: E(eV)=',EEV,' dcxeps=',dcxeps
c***

c add resonances
c energies in eV

         ERES(1)=0.0115
         ERES(2)=0.0315
         ERES(3)=0.0425
         ERES(4)=0.068
         ERES(5)=0.35
         ERES(6)=0.6
         ERES(7)=1.9
         ERES(8)=3.2
         ERES(9)=4.
         ERES(10)=5.
         ERES(11)=7.
         ERES(12)=13.
         ERES(13)=17.
         GAMM(1)=1.226
         GAMM(2)=0.448
         GAMM(3)=0.0965
         GAMM(4)=0.1162
         GAMM(5)=0.6
         GAMM(6)=1.0
         GAMM(7)=0.35
         GAMM(8)=0.25
         GAMM(9)=0.25
         GAMM(10)=0.42
         GAMM(11)=0.28
         GAMM(12)=0.47
         GAMM(13)=0.4
         SRES(1)=54.89
         SRES(2)=7.508
         SRES(3)=8.467
         SRES(4)=4.195
         SRES(5)=2.0
         SRES(6)=8.5
         SRES(7)=0.3
         SRES(8)=0.9
         SRES(9)=0.45
         SRES(10)=1.6
         SRES(11)=0.21
         SRES(12)=0.55
         SRES(13)=0.6
         NRES=13
         DO JR=1,NRES
            X=EEV/ERES(JR)
            DCXEPS=DCXEPS+SRES(JR)/(1.-X**2-CXI*GAMM(JR)*X)
         ENDDO
         DEPS1=REAL(DCXEPS)
         DEPS2=IMAG(DCXEPS)
c*** diagnostic
c         write(0,fmt='(a,1pe10.3,a,1p2e10.3)')
c     &   'free_v4 ckpt 32.2: E(eV)=',EEV,' dcxeps=',dcxeps
c***
      ELSE
         WRITE(0,*)'Fatal error in free: unprepared for IDIEL=',IDIEL
         STOP
      ENDIF
      END
c----------- end of subroutine FREE -------------------------------------
      SUBROUTINE GEOMOPT(REFREL,RAD,WAVE,QABS,QSCA,G)
      IMPLICIT NONE

c Arguments:

      REAL G,QABS,QSCA,RAD,WAVE
      COMPLEX REFREL

c Common variables:

      DOUBLE PRECISION ALPHA,EPS1,EPS22,EPSI1,EPSI22,HALFPI,
     &                 RFRL1,RFRL2,RFRLI1,RFRLI2
      COMMON/GEOMOP/ALPHA,RFRL1,RFRL2,EPS1,EPS22,RFRLI1,RFRLI2,
     &              EPSI1,EPSI22,HALFPI

c Local variables:

      DOUBLE PRECISION DQABS,DQSCAG,FGABS,FGSCAG,QEXT,QRG,X,Z2
      EXTERNAL FGABS,FGSCAG
C***********************************************************************
C Given:
C    RAD = radius of sphere
C    WAVE = wavelength in external medium
C    REFREL = internal (complex) refr. index/external (real) refr. index
C Returns:
C    QABS = absorption cross section/pi*RAD**2 for sphere
C    QSCA = scattering cross section/pi*RAD**2 for sphere
C    G    = <cos(theta)> for scattered radiation
C
C Calculations are made using geometric optics approximations.
C It is assumed that x=2*pi*a/lambda >> 1
C but |m-1|x need not necessarily be large (i.e., not necessarily in
C proper geometric optics regime).  When |m-1|x is not large we make
C simple approximations designed to smoothly bridge the gap between
C Rayleigh-Gans theory for |m-1|x << 1 and proper geometric optics when
C |m-1|x >> 1
C B.T.Draine, Princeton Univ. Observatory, 91.06.17
C History:
C 92.01.06 (BTD) Corrected errors in quadrature over impact parameter.
C 92.01.08 (BTD) Further modifications to smoothly make transition from
C                Rayleigh-Gans to geometric optics.
C                Also added term for 1 internal reflection in FGSCAG
C 92.01.16 (BTD) Changed formula for estimating QEXT (and therefore
C                QSCA=QEXT-QABS) to QEXT=2*QRG/(1.+QRG**2)^{1/2}
C 92.02.24 (BTD) Changed formula for estimating QEXT (and therefore
C                QSCA=QEXT-QABS) to QEXT=QRG/(1.+.25*QRG**2)^{1/2}
C                so that it will be well behaved for QRG << 1 !
C 99.11.01 (BTD) groomed to please FTNCHEK
C end history
C***********************************************************************
      HALFPI=2.D0*ATAN(1.D0)
C
C Compute quantities to go into COMMON/GEOMOP/:
C Quantities required for external reflection:
C
      RFRL1=DBLE(REFREL)
      RFRL2=DBLE(ABS(AIMAG(REFREL)))
      EPS1=RFRL1*RFRL1-RFRL2*RFRL2
      EPS22=(2.D0*RFRL1*RFRL2)**2
C
C Quantities required for internal reflection:
C
      RFRLI1=DBLE(1.E0/REFREL)
      RFRLI2=DBLE(ABS(AIMAG(1.E0/REFREL)))
      EPSI1=RFRLI1*RFRLI1-RFRLI2*RFRLI2
      EPSI22=(2.D0*RFRLI1*RFRLI2)**2
C
C ALPHA=2.*attenuation coefficient*radius
C
      ALPHA=16.D0*HALFPI*RFRL2*DBLE(RAD/WAVE)
      X=4.D0*HALFPI*DBLE(RAD/WAVE)
      Z2=(X*DBLE(ABS(REFREL-1.E0)))**2
C
C It is assumed that this routine will only be invoked if x >> 1.
C It is possible, however, that |m-1|x may not be large.
C We must, therefore, try to produce reasonable behavior when |m-1|x
C is not large.
C
C change 92.01.16:
C
C      QRG=(8./3.)*X*RFRL2+2.*Z2
C
      QRG=(8.D0/3.D0)*X*RFRL2+(32.D0*Z2*X*X/(27.D0+16.D0*X*X))
C
C change 92.02.24:
C      QEXT=2.*QRG/SQRT(1.+QRG**2)
C
      QEXT=QRG/SQRT(1.D0+.25D0*QRG**2)
C
C Note: when (finite) absorption is present, total reflection does not
C occur (cf. Stratton 1941, pp. 500ff if you don't believe it!)
C Perform 16 pt gaussian quadrature
C   Integration variable = [sin(theta)]**2 where theta=angle of incidence
C   of ray relative to surface normal.  Impact parameter = a*sin(theta)
      QABS=0.E0
      IF(ALPHA.GT.0.D0)THEN
          CALL DQG16(0.D0,1.D0,FGABS,DQABS)
          QABS=REAL(DQABS)
      ENDIF
C Make sure that QEXT is at least as large as QABS !
      IF(QEXT.LT.QABS)QEXT=DBLE(QABS)
      CALL DQG16(0.D0,1.D0,FGSCAG,DQSCAG)
C DQSCAG is reflected plus refracted light*cos(theta)
      G=REAL(DQSCAG)/(1.E0-QABS)
C Now need to allow for presence of diffracted light when Q_{ext}>1
      IF(QEXT.GT.1.)THEN
          G=REAL((DQSCAG+QEXT-1.D0)/(QEXT-DBLE(QABS)))
      ENDIF
C Now need to arrive at value for QSCA
      QSCA=REAL(QEXT)-QABS
      RETURN
      END

      FUNCTION FGABS(SINTH2)
      IMPLICIT NONE
C
C Arguments:
C
      DOUBLE PRECISION FGABS,SINTH2
C
C Common:
C
      DOUBLE PRECISION ALPHA,EPS1,EPS22,EPSI1,EPSI22,HALFPI,
     &                 RFRL1,RFRL2,RFRLI1,RFRLI2
      COMMON/GEOMOP/ALPHA,RFRL1,RFRL2,EPS1,EPS22,RFRLI1,RFRLI2,
     &              EPSI1,EPSI22,HALFPI
C
C Local variables:
C
      DOUBLE PRECISION ABSPAR,ABSPER,ATTEN,COSPS,COSTH,NEFF2,
     &                 P2,Q,Q2,R0PAR,R0PER,R1PAR,R1PER,
     &                 SINPS2,TERM1,TERM2,TERM3
C***********************************************************************
C Argument:
C     SINTH2=[sin(theta)]**2 (theta=angle between inc. ray and normal)
C Information in COMMON/GEOMOP/:
C     ALPHA=2*radius*absorption coefficient
C     HALFPI=0.5*pi
C     RFRL1=Re(refractive index of target/medium)
C     RFRL2=Im(refractive index of target/medium)
C     EPS1=Re(epsilon of target/medium)
C     EPS22=[Im(epsilon of target/medium)]**2
C     RFRLI1=Re(refractive index of medium/target)
C     RFRLI2=Im(refractive index of medium/target)
C     EPSI1=Re(epsilon of medium/target)
C     EPSI22=[Im(epsilon of medium/target)]**2
C Returns:
C     FGABS=fraction of incident energy absorbed
C           [where absorbed fraction isaveraged over two
C            incident polarization states]
C It is assumed that material is absorptive (RFRL2>0)
C Expressions taken from Stratton (1941) for nonmagnetic media:
C     NEFF2=[effective refractive index]**2 for Snell's law for
C           refraction from medium into target
C     R0PER=reflection coefficient for perpendicular polarization
C     R0PAR=reflection coefficient for parallel pol.
C     R1PER=internal reflection coefficient for perpendicular pol.
C     R1PAR=internal reflection coefficient for parallel pol.
C
C B.T.Draine, Princeton Univ. Obs., 91.07.17
C history:
C 92.01.06 (BTD) Revised to use different argument (change in variable
C                of integration in calling routine)
C 99.11.01 (BTD) groomed to please FTNCHEK
C end history
C***********************************************************************
C Grazing incidence: assume 100% reflection
      IF(SINTH2.GE.1.D0)THEN
          FGABS=0.D0
          RETURN
      ENDIF
C From here on can assume that cos(theta) is not zero.
C Compute preliminary quantities Q2 and P2
      COSTH=SQRT(1.D0-SINTH2)
      TERM1=EPS1-SINTH2
      TERM2=SQRT(EPS22+TERM1*TERM1)
      Q2=0.5D0*(TERM2+TERM1)
      Q=SQRT(Q2)
      P2=0.5D0*(TERM2-TERM1)
C NEFF2=[effective refractive index]**2 for Snell's law
      NEFF2=Q2+SINTH2
C SINPS2=sin(psi)**2, where psi=angle between refracted ray and normal
      SINPS2=SINTH2/NEFF2
      COSPS=SQRT(1.D0-SINPS2)
C R0PER=reflection coefficient for incident ray, perp. pol.
      R0PER=((Q-COSTH)**2+P2)/((Q+COSTH)**2+P2)
C E parallel to plane defined by incident ray and normal:
C R0PAR=reflection coefficient for incident ray, parallel pol.
      TERM3=SINTH2/COSTH
      R0PAR=R0PER*((Q-TERM3)**2+P2)/((Q+TERM3)**2+P2)
C R1PER=internal reflection coefficient, perpendicular pol.
C R1PAR=internal reflection coefficient, parallel pol.
      TERM1=EPSI1-SINPS2
      TERM2=SQRT(EPSI22+TERM1*TERM1)
      Q2=0.5D0*(TERM2+TERM1)
      Q=SQRT(Q2)
      P2=0.5D0*(TERM2-TERM1)
      R1PER=((Q-COSPS)**2+P2)/((Q+COSPS)**2+P2)
      TERM3=SINPS2/COSPS
      R1PAR=R1PER*((Q-TERM3)**2+P2)/((Q+TERM3)**2+P2)
C ATTEN=attenuation along chord of length 2*RAD*cos(PSI)
      ATTEN=EXP(-ALPHA*COSPS)
C ABSPAR=fraction of "parallel" energy absorbed
C ABSPER=fraction of "perpendicular" energy absorbed
      ABSPAR=(1.D0-R0PAR)*(1.D0-ATTEN)/(1.D0-R1PAR*ATTEN)
      ABSPER=(1.D0-R0PER)*(1.D0-ATTEN)/(1.D0-R1PER*ATTEN)
      FGABS=0.5D0*(ABSPAR+ABSPER)
      RETURN
      END

      FUNCTION FGSCAG(SINTH2)
      IMPLICIT NONE
C Arguments:
      DOUBLE PRECISION FGSCAG,SINTH2
C Common:
      DOUBLE PRECISION ALPHA,EPS1,EPS22,EPSI1,EPSI22,HALFPI,
     &                 RFRL1,RFRL2,RFRLI1,RFRLI2
      COMMON/GEOMOP/ALPHA,RFRL1,RFRL2,EPS1,EPS22,RFRLI1,RFRLI2,
     &              EPSI1,EPSI22,HALFPI
C Local variables:
      DOUBLE PRECISION ATTEN,COSPS,COSTH,NEFF2,P2,PSI,
     &       Q,Q2,R0PAR,R0PER,R1PAR,R1PER,SINPS2,SINTH,
     &       TERM1,TERM2,TERM3,THETA,TPAR1,TPAR2,TPER1,TPER2
C***********************************************************************
C Arguments:
C     SINTH2=[sin(theta)]**2, where theta=angle between incident ray
C            and surface normal
C Information in COMMON/GEOMOP/:
C     ALPHA=2*radius*absorption coefficient
C     HALFPI=0.5*pi
C     RFRL1=Re(refractive index)
C     RFRL2=Im(refractive index)
C     EPS1=RFRL1*RFRL1-RFRL2*RFRL2=Re(epsilon)
C     EPS22=(2*RFRL1*RFRL2)**2=[Im(epsilon)]**2
C Returns:
C     FGSCAG=(scattered fraction)*<cos(theta_s)>
C        where scattered fraction includes direct reflection
C        plus refraction with zero internal reflections,
C        averaged over both incident polarization states.
C
C It is assumed that material is absorptive (RFRL2>0)
C Expressions taken from Stratton (1941) for nonmagnetic media:
C     NEFF2=[effective refractive index]**2 for Snell's law
C     R0PER=reflection coefficient for perpendicular pol.
C     R0PAR=reflection coefficient for parallel pol.
C     R1PER=internal reflection coefficient for perp. pol.
C     R1PAR=internal reflection coefficient for parallel pol.
C Optics is treated in geometrical optics approximation: we assume that
C we may do "ray tracing".
C We assume that only significant contribution to <cos(theta)>
C is due to light either reflected on first incidence or light which
C is refracted once and then propagates without internal reflection.
C
C B.T.Draine, Princeton Univ. Obs., 91.07.17
C history:
C 92.01.06 (BTD) revised: change of argument and use R1PAR,R1PER for
C                internal reflection coefficients
C 99.11.01 (BTD) groomed to please FTNCHEK
C end history
C***********************************************************************
C Grazing incidence: assume 100% reflection
C
      IF(SINTH2.GT.1.D0)THEN
          FGSCAG=1.D0
          RETURN
      ENDIF
C
C From here on can assume that cos(theta) is not zero.
C Compute preliminary quantities Q2 and P2
C
      SINTH=SQRT(SINTH2)
      COSTH=SQRT(1.D0-SINTH2)
      THETA=ASIN(SINTH)
      TERM1=EPS1-SINTH2
      TERM2=SQRT(EPS22+TERM1*TERM1)
      Q2=0.5D0*(TERM2+TERM1)
      Q=SQRT(Q2)
      P2=0.5D0*(TERM2-TERM1)
C
C NEFF2=[effective refractive index]**2 for Snell's law
C
      NEFF2=Q2+SINTH2
C
C SINPSI2=sin(psi)**2, where psi=angle between refracted ray and normal
C
      SINPS2=SINTH2/NEFF2
      COSPS=SQRT(1.D0-SINPS2)
      PSI=ACOS(COSPS)
C
C R0PER=reflection coefficient for incident ray, perp. pol.
C
      R0PER=((Q-COSTH)**2+P2)/((Q+COSTH)**2+P2)
C
C E parallel to plane defined by incident ray and normal:
C R0PAR=reflection coefficient for incdient ray, parallel pol.
C
      TERM3=SINTH2/COSTH
      R0PAR=R0PER*((Q-TERM3)**2+P2)/((Q+TERM3)**2+P2)
C
C R1PER=internal reflection coefficient, perp. pol.
C R1PAR=internal reflection coefficient, parallel pol.
C      TERM1=EPSI1-SINPS2
      TERM2=SQRT(EPSI22+TERM1*TERM1)
      Q2=0.5D0*(TERM2+TERM1)
      Q=SQRT(Q2)
      P2=0.5D0*(TERM2-TERM1)
      R1PER=((Q-COSPS)**2+P2)/((Q+COSPS)**2+P2)
      TERM3=SINPS2/COSPS
      R1PAR=R1PER*((Q-TERM3)**2+P2)/((Q+TERM3)**2+P2)
C
C ATTEN=attenuation along chord of length 2*RAD*cos(PSI)
C
      ATTEN=EXP(-ALPHA*COSPS)
C
C (transmission in)*(transmission out)
C
      TPAR1=(1.D0-R0PAR)*(1.D0-R1PAR)
      TPER1=(1.D0-R0PER)*(1.D0-R1PER)
C
C (transmission in)*(internal reflection)*(transmission out)
C
      TPAR2=TPAR1*R1PAR
      TPER2=TPER1*R1PER
C
C Initial reflection:
C
      FGSCAG=.5D0*(R0PAR+R0PER)*COS(2.D0*(HALFPI-THETA))
C
C Refraction without internal reflection:
C
      FGSCAG=FGSCAG+.5D0*(TPAR1+TPER1)*ATTEN*COS(2.D0*(THETA-PSI))
C
C Refraction plus 1 internal reflection:
C
      FGSCAG=FGSCAG-.5D0*(TPAR2+TPER2)*ATTEN**2*COS(2.D0*THETA-4.D0*PSI)
C
C We assume that higher order terms (two or more internal reflections)
C are negligible because
C (1) very little additional light because of absorption and finite
C     probability of internal reflection
C and/or
C (2) contribution is approximately isotropic.
C
      RETURN
      END
*DECK I1MACH
      INTEGER FUNCTION I1MACH (I)
C***BEGIN PROLOGUE  I1MACH
C***PURPOSE  Return integer machine dependent constants.
C***LIBRARY   SLATEC
C***CATEGORY  R1
C***TYPE      INTEGER (I1MACH-I)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Fox, P. A., (Bell Labs)
C           Hall, A. D., (Bell Labs)
C           Schryer, N. L., (Bell Labs)
C***DESCRIPTION
C
C   I1MACH can be used to obtain machine-dependent parameters for the
C   local machine environment.  It is a function subprogram with one
C   (input) argument and can be referenced as follows:
C
C        K = I1MACH(I)
C
C   where I=1,...,16.  The (output) value of K above is determined by
C   the (input) value of I.  The results for various values of I are
C   discussed below.
C
C   I/O unit numbers:
C     I1MACH( 1) = the standard input unit.
C     I1MACH( 2) = the standard output unit.
C     I1MACH( 3) = the standard punch unit.
C     I1MACH( 4) = the standard error message unit.
C
C   Words:
C     I1MACH( 5) = the number of bits per integer storage unit.
C     I1MACH( 6) = the number of characters per integer storage unit.
C
C   Integers:
C     assume integers are represented in the S-digit, base-A form
C
C                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C                where 0 .LE. X(I) .LT. A for I=0,...,S-1.
C     I1MACH( 7) = A, the base.
C     I1MACH( 8) = S, the number of base-A digits.
C     I1MACH( 9) = A**S - 1, the largest magnitude.
C
C   Floating-Point Numbers:
C     Assume floating-point numbers are represented in the T-digit,
C     base-B form
C                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C                where 0 .LE. X(I) .LT. B for I=1,...,T,
C                0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
C     I1MACH(10) = B, the base.
C
C   Single-Precision:
C     I1MACH(11) = T, the number of base-B digits.
C     I1MACH(12) = EMIN, the smallest exponent E.
C     I1MACH(13) = EMAX, the largest exponent E.
C
C   Double-Precision:
C     I1MACH(14) = T, the number of base-B digits.
C     I1MACH(15) = EMIN, the smallest exponent E.
C     I1MACH(16) = EMAX, the largest exponent E.
C
C   To alter this function for a particular environment, the desired
C   set of DATA statements should be activated by removing the C from
C   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be
C   checked for consistency with the local operating system.
C
C***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
C                 a portable library, ACM Transactions on Mathematical
C                 Software 4, 2 (June 1978), pp. 177-188.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   750101  DATE WRITTEN
C   891012  Added VAX G-floating constants.  (WRB)
C   891012  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900618  Added DEC RISC constants.  (WRB)
C   900723  Added IBM RS 6000 constants.  (WRB)
C   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16.
C           (RWC)
C   910710  Added HP 730 constants.  (SMR)
C   911114  Added Convex IEEE constants.  (WRB)
C   920121  Added SUN -r8 compiler option constants.  (WRB)
C   920229  Added Touchstone Delta i860 constants.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   920625  Added Convex -p8 and -pd8 compiler option constants.
C           (BKS, WRB)
C   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
C   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler
C           options.  (DWL, RWC and WRB).
C***END PROLOGUE  I1MACH
C
      INTEGER IMACH(16),OUTPUT
      SAVE IMACH
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT COMPILER
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1022 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        129 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1025 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA IMACH( 1) /          7 /
C     DATA IMACH( 2) /          2 /
C     DATA IMACH( 3) /          2 /
C     DATA IMACH( 4) /          2 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         33 /
C     DATA IMACH( 9) / Z1FFFFFFFF /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -256 /
C     DATA IMACH(13) /        255 /
C     DATA IMACH(14) /         60 /
C     DATA IMACH(15) /       -256 /
C     DATA IMACH(16) /        255 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         48 /
C     DATA IMACH( 6) /          6 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /          8 /
C     DATA IMACH(11) /         13 /
C     DATA IMACH(12) /        -50 /
C     DATA IMACH(13) /         76 /
C     DATA IMACH(14) /         26 /
C     DATA IMACH(15) /        -50 /
C     DATA IMACH(16) /         76 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         48 /
C     DATA IMACH( 6) /          6 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /          8 /
C     DATA IMACH(11) /         13 /
C     DATA IMACH(12) /        -50 /
C     DATA IMACH(13) /         76 /
C     DATA IMACH(14) /         26 /
C     DATA IMACH(15) /     -32754 /
C     DATA IMACH(16) /      32780 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          8 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 9223372036854775807 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /      -4095 /
C     DATA IMACH(13) /       4094 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /      -4095 /
C     DATA IMACH(16) /       4094 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /    6LOUTPUT/
C     DATA IMACH( 5) /         60 /
C     DATA IMACH( 6) /         10 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         48 /
C     DATA IMACH( 9) / 00007777777777777777B /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /       -929 /
C     DATA IMACH(13) /       1070 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /       -929 /
C     DATA IMACH(16) /       1069 /
C
C     MACHINE CONSTANTS FOR THE CELERITY C1260
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          0 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / Z'7FFFFFFF' /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1022 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fn COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fi COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -p8 COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 9223372036854775807 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         53 /
C     DATA IMACH(12) /      -1023 /
C     DATA IMACH(13) /       1023 /
C     DATA IMACH(14) /        113 /
C     DATA IMACH(15) /     -16383 /
C     DATA IMACH(16) /      16383 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -pd8 COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 9223372036854775807 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         53 /
C     DATA IMACH(12) /      -1023 /
C     DATA IMACH(13) /       1023 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE CRAY
C     USING THE 46 BIT INTEGER COMPILER OPTION
C
C     DATA IMACH( 1) /        100 /
C     DATA IMACH( 2) /        101 /
C     DATA IMACH( 3) /        102 /
C     DATA IMACH( 4) /        101 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          8 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         46 /
C     DATA IMACH( 9) / 1777777777777777B /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /      -8189 /
C     DATA IMACH(13) /       8190 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /      -8099 /
C     DATA IMACH(16) /       8190 /
C
C     MACHINE CONSTANTS FOR THE CRAY
C     USING THE 64 BIT INTEGER COMPILER OPTION
C
C     DATA IMACH( 1) /        100 /
C     DATA IMACH( 2) /        101 /
C     DATA IMACH( 3) /        102 /
C     DATA IMACH( 4) /        101 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          8 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 777777777777777777777B /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /      -8189 /
C     DATA IMACH(13) /       8190 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /      -8099 /
C     DATA IMACH(16) /       8190 /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     DATA IMACH( 1) /         11 /
C     DATA IMACH( 2) /         12 /
C     DATA IMACH( 3) /          8 /
C     DATA IMACH( 4) /         10 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /         16 /
C     DATA IMACH(11) /          6 /
C     DATA IMACH(12) /        -64 /
C     DATA IMACH(13) /         63 /
C     DATA IMACH(14) /         14 /
C     DATA IMACH(15) /        -64 /
C     DATA IMACH(16) /         63 /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING G_FLOAT
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING IEEE_FLOAT
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE DEC RISC
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     USING D_FLOATING
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     USING G_FLOATING
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         32 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1022 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          0 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         24 /
C     DATA IMACH( 6) /          3 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         23 /
C     DATA IMACH( 9) /    8388607 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         23 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         38 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /         43 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          6 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         63 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HP 730
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          4 /
C     DATA IMACH( 4) /          1 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         23 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         39 /
C     DATA IMACH(15) /       -128 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          4 /
C     DATA IMACH( 4) /          1 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         23 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         55 /
C     DATA IMACH(15) /       -128 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          7 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         32 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1015 /
C     DATA IMACH(16) /       1017 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) /  Z7FFFFFFF /
C     DATA IMACH(10) /         16 /
C     DATA IMACH(11) /          6 /
C     DATA IMACH(12) /        -64 /
C     DATA IMACH(13) /         63 /
C     DATA IMACH(14) /         14 /
C     DATA IMACH(15) /        -64 /
C     DATA IMACH(16) /         63 /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C
      DATA IMACH( 1) /          5 /
      DATA IMACH( 2) /          6 /
      DATA IMACH( 3) /          0 /
      DATA IMACH( 4) /          0 /
      DATA IMACH( 5) /         32 /
      DATA IMACH( 6) /          4 /
      DATA IMACH( 7) /          2 /
      DATA IMACH( 8) /         31 /
      DATA IMACH( 9) / 2147483647 /
      DATA IMACH(10) /          2 /
      DATA IMACH(11) /         24 /
      DATA IMACH(12) /       -125 /
      DATA IMACH(13) /        127 /
      DATA IMACH(14) /         53 /
      DATA IMACH(15) /      -1021 /
      DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE IBM RS 6000
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          0 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE INTEL i860
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          5 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         54 /
C     DATA IMACH(15) /       -101 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          5 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         62 /
C     DATA IMACH(15) /       -128 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE SUN
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE SUN
C     USING THE -r8 COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         53 /
C     DATA IMACH(12) /      -1021 /
C     DATA IMACH(13) /       1024 /
C     DATA IMACH(14) /        113 /
C     DATA IMACH(15) /     -16381 /
C     DATA IMACH(16) /      16384 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          1 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         60 /
C     DATA IMACH(15) /      -1024 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
C
C     DATA IMACH( 1) /          1 /
C     DATA IMACH( 2) /          1 /
C     DATA IMACH( 3) /          0 /
C     DATA IMACH( 4) /          1 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C***FIRST EXECUTABLE STATEMENT  I1MACH
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10
C
      I1MACH = IMACH(I)
      RETURN
C
   10 CONTINUE
      WRITE (UNIT = OUTPUT, FMT = 9000)
 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
C
C     CALL FDUMP
C
      STOP
      END
      SUBROUTINE INDEX(A,TGR,WAVE,IDIEL,EPS11,EPS2,ENRE1,ENIM)
      IMPLICIT NONE
c------------------------------ index_v6 -------------------------------
C Arguments:

      INTEGER IDIEL
      REAL A,ENIM,ENRE1,EPS11,EPS2,TGR,WAVE

C Parameters:

      INTEGER NTABMX,NTYPES
      PARAMETER(NTABMX=10000,NTYPES=50)

c-----------------------------------------------------------------------

c Common variables:

      INTEGER NTAB
      REAL EN,ENREAL,ENIMAG
      COMMON/TABCOM2/EN(1:NTABMX),ENREAL(1:NTABMX),ENIMAG(1:NTABMX),NTAB
c-----------------------------------------------------------------------
      REAL APORO,FPORO
      COMMON/POROSITY/APORO,FPORO
c-----------------------------------------------------------------------

C Local variables:

C Note: length of PATH must agree with specification in FINDPATH

      CHARACTER*15 PATH

      INTEGER IDIELL,IDIELT,ITHEORY,J
      INTEGER NTABST(NTYPES)
      REAL ABSCOF,ABSIND,DDEPS1,DDEPS2,ENIM0,ENRE,ENRE10,EX,X,Y,Z
      REAL ENST(NTABMX,NTYPES),ENREST(NTABMX,NTYPES),
     &     ENIMST(NTABMX,NTYPES)
      DOUBLE PRECISION DEPS1,DEPS2,DROOT,PORO,XPORO
      DOUBLE COMPLEX CXEPS,CXEPSEFF,CXEPSINC
      CHARACTER*40 FILEIN(NTYPES)

      INTRINSIC DIMAG

      EXTERNAL FREE

      DATA APORO/0.15/,FPORO/0.50/
C***********************************************************************
C
C Subroutine to provide dielectric constants of selected materials
C Input:
C   A=grain radius (micron)
C   TGR=grain temperature (degK)
C   WAVE=wavelength (micron)
C   IDIEL=1 = DL84 graphite, E || c , sigma_dc=29.1 mho/cm
C         2 = DL84 graphite, E perpendicular to c axis
C         3 = astronomical silicate from DL84
C         4 = Type Ib diamond with N/C=0.003
C         5 = meteoritic diamond (Allende CJ)
C         6 = alpha-SiC
C         7 = "smoothed UV" 2000 astronomical silicate
C         8 = "new" 2000 astronomical silicate
C         9 = graphite, E || c axis, with epsilon(DL84) except for
C               increased DC conductivity (same as IDIEL=1 but with
C               sigma_dc=62.5 mho/cm instead of 29.1 mho/cm)
C         10 = crystalline H2O ice
C         11 = D03 astronomical silicate
C         12 = D03 graphite E || c ,
c                with enhanced conductivity sigma_dc = 62.5 mho/cm
C         13 = D03 graphite E perp to c
C         14 = D03 astronomical silicate with porosity(a)
C         15 = D03 graphite, E parallel to c, with porosity(a)
c                with enhanced DC conductivity
C         16 = D03 graphite, E perp. to c, with porosity(a)
C         17 = Cpaexp1: E parallel to c. Djuric & Li 99, optical data
C         18 = Cpaexp2: E parallel to c. Djuric & Li 99, EELS data
C         19 = Cpeexp : E perp to c. Djuric & Li 99
C         20 = 800C pyrolyzed cellulose (Jaeger etal 98) extrap to submm
C         21 = 600C pyrolyzed cellulose (Jaeger etal 98) extrap to submm
C         22 = amorph.car. AC1 from Rouleau & Martin 1991
C         23 = amorph.car. BE1 from Rouleau & Martin 1991
C         24 = amorph.car. FC21PS from Rouleau & Martin 1991
C         25 = amorph.car. HAPS from Rouleau & Martin 1991
C         26 = amorph.car. BE from Zubko etal 1996
C         27 = amorph.car. ACAR from Zubko etal 1996
C         28 = amorph.car. ACH2 from Zubko etal 1996
C         29 = amorph.car amcACH2 based on ACH2 from Zubko etal 1996
C                           but revised above 7 eV
C         30 = D03 graphite, E || c, but with a two-component model for
C                           the free electron contribution, to suppress
C                           the pronounced peak in Q_abs at 33um
c                           which occurred for the DL07 parameters
c         31 = Fe metal
c         32 = magnetite Fe3O4
c         33 = maghemite gamma-Fe2O3
c
c Note: for IDIEL=14,15,16, porosity parameters APORO,FPORO are set at
c       compile time by DATA statement 
c       However, user can reset APORO,FPORO [if desired]
c       by passing them through COMMON/POROSITY/

C Returns:
C   EPS11 = Re(epsilon)-1 (epsilon=complex dielectric constant)
C   EPS2  = Im(epsilon)
C   ENRE1 = Re(m)-1  (m=complex refractive index)
C   ENIM  = Im(m)
C
C Uses following subroutines:
C   FREE to compute free electron contribtion to dielectric function
C       (only point where temperature/size dependence is now present)
C   TABLE to interpolate in table of optical constants
C       PARAB3,PARABA (called by TABLE)
C       PARAB4,PARABB (called by TABLE)
C   REFICE to provide optical constants of H2O ice.
C
c Requires following data files:
c  1   '/u/draine/work/opt/dat/index.Cpa'
c  2   '/u/draine/work/opt/dat/index.Cpe'
c  3   '/u/draine/work/opt/dat/index.sil'
c  4   '/u/draine/work/opt/dat/index.diamIb'
c  5   '/u/draine/work/opt/dat/index.mdiam'
c  6   '/u/draine/work/opt/dat/index.SiC'
c  7   '/u/draine/work/opt/dat/index.sil_suv'
c  8   '/u/draine/work/opt/dat/index.sil_new'
c  9   same as 1
c 10   none
c 11   '/u/draine/work/opt/dat/index_silD03'
c 12   '/u/draine/work/opt/dat/index_CpaD03'
c 13   '/u/draine/work/opt/dat/index_CpeD03'
c 14   same as 11
c 15   same as 12
c 16   same as 13
c 17   '/u/draine/work/opt/dat/index_Cpaexp1'
c 18   '/u/draine/work/opt/dat/index_Cpaexp2'
c 19   '/u/draine/work/opt/dat/index_Cpaexp'
c 20   '/u/draine/work/opt/dat/index_pcel800'
c 21   '/u/draine/work/opt/dat/index_pcel600'
c 22   '/u/draine/work/opt/dat/index_amcAC1_RM91'
c 23   '/u/draine/work/opt/dat/index_amcBE1_RM91'
c 24   '/u/draine/work/opt/dat/index_amcFC21PS_RM91'
c 25   '/u/draine/work/opt/dat/index_amcHAPS_RM91'
c 26   '/u/draine/work/opt/dat/index_amcBE_ZMCB96'
c 27   '/u/draine/work/opt/dat/index_amcACAR_ZMCB96'
c 28   '/u/draine/work/opt/dat/index_amcACH2_ZMCB96'
c 29   '/u/draine/work/opt/dat/index_amcACH2'
c 30   same as 12
c 31   '/u/draine/work/opt/dat/index_Fe'
c 32   '/u/draine/work/opt/dat/index_fe3o4'
c 33   '/u/draine/work/opt/dat/index_fe2o3'
c
c for IDIEL=14,15,16 porosity is determined by parameters APORO,FPORO
c which are specified in DATA statement above
c values of APORO,FPORO can be read out via COMMON/POROSITY/
c and value can be changed via COMMON/POROSITY/
c 
c for a < APORO, void fraction = 0 (compact grain)
c     a > APORO, void fraction = FPORO*[1-(APORO/a)^2]
c                             --> FPORO in limit a >> APORO

C***********************************************************************
C B.T. Draine, Princeton Univ. Obs.,
C History:
C 89.01.20 (BTD): Modified...
C 89.02.01 (BTD): Modified to do double precision arithmetic in
C                 computing ENRE,ENIM
C 91.07.21 (BTD): Explicit declaration of all variables.
C 91.12.12 (BTD): Added SiC to materials included.
C 91.12.23 (BTD): Modified to return Re(m-1) and Im(m) rather than
C                 Re(m) and Im(m), in order to deal with materials with
C                 |m-1| << 1 (e.g., hard X-ray region)
C 99.02.17 (BTD): Add SAVE ENIMST,ENREST,ENST,IDIELL,NTABST to
C                 make code f77 compliant (and compatible with g77)
C 99.02.18 (BTD): Modified to allow new option of "smoothed UV silicate"
C 00.06.27 (BTD): Put trap in to halt if called with IDIEL=-1 or >10
C 00.07.17 (BTD): Modified to support IDIEL=8
C 00.09.13 (BTD): Modified to support IDIEL=9 (new graphite, E par.c)
C                 introduced new local variable IDIELT to facilitate this
C 02.05.18 (BTD): Modified to support IDIEL=17 
C                 silD02 with improved treatment of X-ray absorption
C                 changed NTYPE to 20
C 02.05.19 (BTD): Modified to support IDIEL=18 and 19
C                 CpaD02 and CpeD02 graphite with improved treatment of
C                 X-ray absorption
C 02.12.22 (BTD): Modified to support IDIEL=20
C                 silD02 material in a porous grain
C                 with porosity depending on grain radius
C                 parameters APORO and FPORO determine porosity of grain
C                 of radius A
C 03.01.02 (BTD): Changed IDIEL option numbers: 22-24 are now porous
C                 silicate and graphite options.
C 05.09.19 (BTD): Added variable PATH and call to FINDPATH to permit
C                 use with different directory structures
C                 Note that subroutine FINDPATH must set appropriate
C                 path to directory where opt/dat resides
C 07.12.20 (BTD): v3
C                 Add option IDIEL=44 for D03 graphite, E||c, but with
C                 a two-component model for the free electron
C                 contribution to epsilon in order to suppress the
C                 artificial peak at 33um that resulted from the
C                 free electron contribution for the free-electron
C                 parameters used by DL07.
c 07.12.21 (BTD): v4
C                 Renumber dielectric options
c 11.02.13 (BTD): v5
c                 Add option IDIEL=31 for metallic Fe
c 11.09.11 (BTD): v6
c                 Add support for options
c                 IDIEL=32 for magnetite Fe3O4
c                 IDIEL=33 for maghemite gamma-Fe2O3
C end history
C***********************************************************************
      SAVE ENIMST,ENREST,ENST,FILEIN,IDIELL,NTABST
      DATA IDIELL/0/,NTABST/NTYPES*0/
      DATA (FILEIN(J),J=1,8)/
     &'opt/dat/index.Cpa',
     &'opt/dat/index.Cpe',
     &'opt/dat/index.sil',
     &'opt/dat/index.diamIb',
     &'opt/dat/index.mdiam',
     &'opt/dat/index.SiC',
     &'opt/dat/index.sil_suv',
     &'opt/dat/index.sil_new'/

      DATA (FILEIN(J),J=11,13)/
     &'opt/dat/index_silD03',
     &'opt/dat/index_CpaD03',
     &'opt/dat/index_CpeD03'/

      DATA (FILEIN(J),J=17,21)/
     &'opt/dat/index_Cpaexp1',
     &'opt/dat/index_Cpaexp2',
     &'opt/dat/index_Cpeexp',
     &'opt/dat/index_pcel800',
     &'opt/dat/index_pcel600'/

      DATA (FILEIN(J),J=22,25)/
     &'opt/dat/index_amcAC1_RM91',
     &'opt/dat/index_amcBE1_RM91',
     &'opt/dat/index_amcFC21PS_RM91',
     &'opt/dat/index_amcHAPS_RM91'/

      DATA (FILEIN(J),J=26,29)/
     &'opt/dat/index_amcBE_ZMCB96',
     &'opt/dat/index_amcACAR_ZMCB96',
     &'opt/dat/index_amcACH2_ZMCB96',
     &'opt/dat/index_amcACH2'/

      DATA (FILEIN(J),J=31,33)/
     &'opt/dat/index_Fe',
     &'opt/dat/index_fe3o4',
     &'opt/dat/index_fe2o3'/

      IF(IDIEL.EQ.10)THEN
C H2O ice properties are provided by routine REFICE
          CALL REFICE(0,WAVE,TGR,ENRE,ENIM,ABSIND,ABSCOF)
          ENRE1=ENRE-1.
          EPS11=ENRE*ENRE-ENIM*ENIM-1.
          EPS2=2.*ENRE*ENIM
          RETURN
      ENDIF

C Proceed to use tables

      IDIELT=IDIEL

c special cases:

      IF(IDIEL.EQ.9)IDIELT=1   ! DL84 graph(E||c)+increased sig_d + poro
      IF(IDIEL.EQ.14)IDIELT=11 ! D03 silicate + porosity
      IF(IDIEL.EQ.15)IDIELT=12 ! D03 graph(E||c)+increased sig_dc + poro
      IF(IDIEL.EQ.16)IDIELT=13 ! D03 graph(E perp c) + porosity
      IF(IDIEL.EQ.30)IDIELT=12

C***** Check if desired optical constants are already in /TABCOM2/:

      IF(IDIELL.EQ.IDIELT)GOTO 2000

C***** Define arrays EN,ENREAL,ENIMAG in /TABCOM2/:
C Determine if have already read needed file:

      IF(NTABST(IDIELT).GT.0)GOTO 1500

C Have not yet read needed index file: open appropriate file

      CALL FINDPATH(PATH)

      OPEN(UNIT=2,FILE=PATH//FILEIN(IDIELT),STATUS='OLD')

C Read file:
C Skip first two lines in file (comments)

      READ(2,0099)
 0099 FORMAT(/)
      J=0

C Each line of file contains: E(eV), Re(n-1), Im(n)
C (where n = complex refractive index)

 0100 READ(2,*,END=0199)X,Y,Z
      J=J+1
c*** diagnostic
c      write(0,*)'index_v6 ckpt 10: j=',j,' x=',x
c**
      IF(J.GT.NTABMX)THEN
         WRITE(0,*)'Fatal error in index: input file',FILEIN(IDIELT),
     &   ' has more than NTABMX=',NTABMX,' lines'
         STOP
      ENDIF
      ENST(J,IDIELT)=X
      ENREST(J,IDIELT)=Y
      ENIMST(J,IDIELT)=Z
      GOTO 0100
 0199 NTABST(IDIELT)=J
c*** diagnostic
c      write(0,*)'index_v6 ckp 11'
c***
      CLOSE(UNIT=2,STATUS='KEEP')

C Now redefine arrays in /TABCOM2/:

 1500 NTAB=NTABST(IDIELT)
      DO 1700 J=1,NTAB
         EN(J)=ENST(J,IDIELT)
         ENREAL(J)=ENREST(J,IDIELT)

C Beware of differing conventions regarding sign of Im(epsilon): here
C require it to be positive.

         ENIMAG(J)=ABS(ENIMST(J,IDIELT))
 1700 CONTINUE
      CONTINUE
      IDIELL=IDIELT

c*** diagnostic: check data integrity
c      do j=1,ntab
c         if(enimag(j).lt.0.)then
c            write(0,*)' problem in index, loading tabcom2:'
c            write(0,*)' j=',j,' en(j)=',en(j),' enimag(j)=',enimag(j)
c         endif
c      enddo
c***

C*** End initialization of /TABCOM2/
C    Subroutine TABLE is now ready for use

 2000 EX=1.23984E0/WAVE

C EX=energy (eV)
C WAVE=wavelength (micron)

c*** diagnostic: check data integrity
c      do j=1,ntab
c         if(enimag(j).lt.0.)then
c            write(0,*)' problem in index: tabcom2 has been corrupted'
c            write(0,*)' j=',j,' en(j)=',en(j),' enimag(j)=',enimag(j)
c         endif
c      enddo
c***
      CALL TABLE2(EX,ENRE10,ENIM0)

c*** diagnostic
c      if(enim0.lt.0.)then
c         write(0,*)'Problem in INDEX: TABLE returned ENIM0=',ENIM0
c         write(0,*)'   for EX   =',EX
c         write(0,*)'       WAVE =',WAVE
c         write(0,*)'       IDIEL=',IDIEL
c      endif
c***

C 1+ENRE10+i*ENIM0=refractive index due to bound electrons
C 1+DEPS1+i*DEPS2=dielectric function

      DEPS1=DBLE(ENRE10*(ENRE10+2.E0)-ENIM0**2)
      DEPS2=2.D0*DBLE(ENIM0*(1.+ENRE10))

C Add contribution of conduction band electrons

c*** diagnostic
c      write(0,*)'about to call free with ex=',ex
c***
      CALL FREE(EX,IDIEL,TGR,A,DDEPS1,DDEPS2)

c*** diagnostic
      if(ddeps2.lt.0.)then
         write(0,*)'index ckpt 3: returned from free with ddeps2=',
     &             ddeps2
      endif
c***
      
      DEPS1=DEPS1+DBLE(DDEPS1)
      DEPS2=DEPS2+DBLE(DDEPS2)

C For IDIEL=22,23,24, assume size-dependent porosity and use effective 
C medium theory to calculate effective dielectric function

      IF(IDIEL.EQ.14.OR.
     &   IDIEL.EQ.15.OR.
     &   IDIEL.EQ.16)THEN

c APORO = radius (um) below which grain is assumed to have zero porosity
c FPORO = porosity in limit of A >> APORO
c         (porosity = filling factor of vacuum)
c ITHEORY = 0 for simple linear average of dielectric functions
c           1 for Bruggemann theory
c           2 for Garnett theory

         APORO=0.2
         FPORO=0.5
         ITHEORY=1
         IF(A.GT.APORO)THEN

c model: when grain is larger than APORO, additional volume is
c        assumed to have porosity FPORO

            XPORO=(1.-(APORO/A)**3)
            PORO=FPORO*XPORO
            CXEPS=1.D0+DEPS1+(0.D0,1.D0)*DEPS2
            CXEPSINC=(1.D0,0.D0)
            CALL EFFMED(ITHEORY,PORO,CXEPS,CXEPSINC,CXEPSEFF)
            DEPS1=DBLE(CXEPSEFF)-1.D0
            DEPS2=DIMAG(CXEPSEFF)
         ENDIF
      ENDIF

C Compute new complex refractive index

      IF((DEPS1**2+DEPS2**2).LT.1.D-6)THEN

C Both DEPS1 and DEPS2 small:
C Expand to 2nd order in powers of DEPS1 and DEPS2:

         ENRE1=REAL(.5D0*DEPS1-.125D0*(DEPS1**2-DEPS2**2))
         ENIM=REAL(.5D0*DEPS2-.25D0*DEPS1*DEPS2)
      ELSEIF(DEPS2.LT.1.D-3*ABS(1.D0+DEPS1))THEN
C
C DEPS1 not small, but DEPS2/DEPS1 small:
C
         ENRE1=REAL(SQRT(1.D0+DEPS1)*
     &         (1.D0+.125D0*(DEPS2/(1.D0+DEPS1))**2)-1.D0)
         ENIM=REAL(SQRT(1.D0+DEPS1)*.5D0*(DEPS2/(1.D0+DEPS1)))
      ELSE
         DROOT=SQRT((1.D0+DEPS1)**2+DEPS2**2)
         ENRE1=REAL(SQRT(.5D0*(DROOT+1.D0+DEPS1))-1.D0)
         ENIM=REAL(SQRT(.5D0*(DROOT-1.D0-DEPS1)))
      ENDIF

C Compute single precision eps1,eps2:

      EPS11=REAL(DEPS1)
      EPS2=REAL(DEPS2)

c*** sanity check
      IF(EPS2.LT.0.)THEN
         WRITE(0,*)'Fatal error in subroutine INDEX:'
         WRITE(0,*)'  EPS2=',EPS2,' for IDIEL=',IDIEL
         WRITE(0,*)'  A=',A,' TGR=',TGR,' WAVE=',WAVE
         STOP
      ENDIF
c***

      RETURN
      END





      DOUBLE PRECISION FUNCTION BESJ0(X)
CS    REAL FUNCTION BESJ0(X)
C--------------------------------------------------------------------
C
C This subprogram computes approximate values for Bessel functions
C   of the first kind of order zero for arguments  |X| <= XMAX
C   (see comments heading CALJY0).
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL  X, RESULT
      DOUBLE PRECISION  X, RESULT
C--------------------------------------------------------------------
      JINT=0
      CALL CALJY0(X,RESULT,JINT)
      BESJ0 = RESULT
      RETURN
C---------- Last line of BESJ0 ----------
      END
      DOUBLE PRECISION FUNCTION BESY0(X)
CS    REAL FUNCTION BESY0(X)
C--------------------------------------------------------------------
C
C This subprogram computes approximate values for Bessel functions
C   of the second kind of order zero for arguments 0 < X <= XMAX
C   (see comments heading CALJY0).
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL  X, RESULT
      DOUBLE PRECISION  X, RESULT
C--------------------------------------------------------------------
      JINT=1
      CALL CALJY0(X,RESULT,JINT)
      BESY0 = RESULT
      RETURN
C---------- Last line of BESY0 ----------
      END
      SUBROUTINE CALJY0(ARG,RESULT,JINT)
C---------------------------------------------------------------------
C
C This packet computes zero-order Bessel functions of the first and
C   second kind (J0 and Y0), for real arguments X, where 0 < X <= XMAX
C   for Y0, and |X| <= XMAX for J0.  It contains two function-type
C   subprograms,  BESJ0  and  BESY0,  and one subroutine-type
C   subprogram,  CALJY0.  The calling statements for the primary
C   entries are:
C
C           Y = BESJ0(X)
C   and
C           Y = BESY0(X),
C
C   where the entry points correspond to the functions J0(X) and Y0(X),
C   respectively.  The routine  CALJY0  is intended for internal packet
C   use only, all computations within the packet being concentrated in
C   this one routine.  The function subprograms invoke  CALJY0  with
C   the statement
C           CALL CALJY0(ARG,RESULT,JINT),
C   where the parameter usage is as follows:
C
C      Function                  Parameters for CALJY0
C       call              ARG             RESULT          JINT
C
C     BESJ0(ARG)     |ARG| .LE. XMAX       J0(ARG)          0
C     BESY0(ARG)   0 .LT. ARG .LE. XMAX    Y0(ARG)          1
C
C   The main computation uses unpublished minimax rational
C   approximations for X .LE. 8.0, and an approximation from the 
C   book  Computer Approximations  by Hart, et. al., Wiley and Sons, 
C   New York, 1968, for arguments larger than 8.0   Part of this
C   transportable packet is patterned after the machine-dependent
C   FUNPACK program BESJ0(X), but cannot match that version for
C   efficiency or accuracy.  This version uses rational functions
C   that are theoretically accurate to at least 18 significant decimal
C   digits for X <= 8, and at least 18 decimal places for X > 8.  The
C   accuracy achieved depends on the arithmetic system, the compiler,
C   the intrinsic functions, and proper selection of the machine-
C   dependent constants.
C
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   XINF   = largest positive machine number
C   XMAX   = largest acceptable argument.  The functions AINT, SIN
C            and COS must perform properly for  ABS(X) .LE. XMAX.
C            We recommend that XMAX be a small integer multiple of
C            sqrt(1/eps), where eps is the smallest positive number
C            such that  1+eps > 1. 
C   XSMALL = positive argument such that  1.0-(X/2)**2 = 1.0
C            to machine precision for all  ABS(X) .LE. XSMALL.
C            We recommend that  XSMALL < sqrt(eps)/beta, where beta
C            is the floating-point radix (usually 2 or 16).
C
C     Approximate values for some important machines are
C
C                          eps      XMAX     XSMALL      XINF  
C
C  CDC 7600      (S.P.)  7.11E-15  1.34E+08  2.98E-08  1.26E+322
C  CRAY-1        (S.P.)  7.11E-15  1.34E+08  2.98E-08  5.45E+2465
C  IBM PC (8087) (S.P.)  5.96E-08  8.19E+03  1.22E-04  3.40E+38
C  IBM PC (8087) (D.P.)  1.11D-16  2.68D+08  3.72D-09  1.79D+308
C  IBM 195       (D.P.)  2.22D-16  6.87D+09  9.09D-13  7.23D+75
C  UNIVAC 1108   (D.P.)  1.73D-18  4.30D+09  2.33D-10  8.98D+307
C  VAX 11/780    (D.P.)  1.39D-17  1.07D+09  9.31D-10  1.70D+38
C
C*******************************************************************
C*******************************************************************
C
C Error Returns
C
C  The program returns the value zero for  X .GT. XMAX, and returns
C    -XINF when BESLY0 is called with a negative or zero argument.
C
C
C Intrinsic functions required are:
C
C     ABS, AINT, COS, LOG, SIN, SQRT
C
C
C  Latest modification: June 2, 1989
C
C  Author: W. J. Cody
C          Mathematics and Computer Science Division 
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C--------------------------------------------------------------------
      INTEGER I,JINT
CS    REAL
      DOUBLE PRECISION
     1       ARG,AX,CONS,DOWN,EIGHT,FIVE5,FOUR,ONE,ONEOV8,PI2,PJ0,
     2       PJ1,PLG,PROD,PY0,PY1,PY2,P0,P1,P17,QJ0,QJ1,QLG,QY0,QY1,
     3       QY2,Q0,Q1,RESJ,RESULT,R0,R1,SIXTY4,THREE,TWOPI,TWOPI1,
     4       TWOPI2,TWO56,UP,W,WSQ,XDEN,XINF,XMAX,XNUM,XSMALL,XJ0,
     5       XJ1,XJ01,XJ02,XJ11,XJ12,XY,XY0,XY01,XY02,XY1,XY11,XY12,
     6       XY2,XY21,XY22,Z,ZERO,ZSQ
      DIMENSION PJ0(7),PJ1(8),PLG(4),PY0(6),PY1(7),PY2(8),P0(6),P1(6),
     1          QJ0(5),QJ1(7),QLG(4),QY0(5),QY1(6),QY2(7),Q0(5),Q1(5)
C-------------------------------------------------------------------
C  Mathematical constants
C    CONS = ln(.5) + Euler's gamma
C-------------------------------------------------------------------
CS    DATA ZERO,ONE,THREE,FOUR,EIGHT/0.0E0,1.0E0,3.0E0,4.0E0,8.0E0/,
CS   1     FIVE5,SIXTY4,ONEOV8,P17/5.5E0,64.0E0,0.125E0,1.716E-1/,
CS   2     TWO56,CONS/256.0E0,-1.1593151565841244881E-1/,
CS   3     PI2,TWOPI/6.3661977236758134308E-1,6.2831853071795864769E0/,
CS   4     TWOPI1,TWOPI2/6.28125E0,1.9353071795864769253E-3/
      DATA ZERO,ONE,THREE,FOUR,EIGHT/0.0D0,1.0D0,3.0D0,4.0D0,8.0D0/,
     1     FIVE5,SIXTY4,ONEOV8,P17/5.5D0,64.0D0,0.125D0,1.716D-1/,
     2     TWO56,CONS/256.0D0,-1.1593151565841244881D-1/,
     3     PI2,TWOPI/6.3661977236758134308D-1,6.2831853071795864769D0/,
     4     TWOPI1,TWOPI2/6.28125D0,1.9353071795864769253D-3/
C-------------------------------------------------------------------
C  Machine-dependent constants
C-------------------------------------------------------------------
CS    DATA XMAX/8.19E+03/,XSMALL/1.22E-09/,XINF/1.7E+38/
      DATA XMAX/1.07D+09/,XSMALL/9.31D-10/,XINF/1.7D+38/
C-------------------------------------------------------------------
C  Zeroes of Bessel functions
C-------------------------------------------------------------------
CS    DATA XJ0/2.4048255576957727686E+0/,XJ1/5.5200781102863106496E+0/,
CS   1     XY0/8.9357696627916752158E-1/,XY1/3.9576784193148578684E+0/,
CS   2     XY2/7.0860510603017726976E+0/,
CS   3     XJ01/ 616.0E+0/, XJ02/-1.4244423042272313784E-03/,
CS   4     XJ11/1413.0E+0/, XJ12/ 5.4686028631064959660E-04/,
CS   5     XY01/ 228.0E+0/, XY02/ 2.9519662791675215849E-03/,
CS   6     XY11/1013.0E+0/, XY12/ 6.4716931485786837568E-04/,
CS   7     XY21/1814.0E+0/, XY22/ 1.1356030177269762362E-04/
      DATA XJ0/2.4048255576957727686D+0/,XJ1/5.5200781102863106496D+0/,
     1     XY0/8.9357696627916752158D-1/,XY1/3.9576784193148578684D+0/,
     2     XY2/7.0860510603017726976D+0/,
     3     XJ01/ 616.0D+0/, XJ02/-1.4244423042272313784D-03/,
     4     XJ11/1413.0D+0/, XJ12/ 5.4686028631064959660D-04/,
     5     XY01/ 228.0D+0/, XY02/ 2.9519662791675215849D-03/,
     6     XY11/1013.0D+0/, XY12/ 6.4716931485786837568D-04/,
     7     XY21/1814.0D+0/, XY22/ 1.1356030177269762362D-04/
C-------------------------------------------------------------------
C  Coefficients for rational approximation to ln(x/a)
C--------------------------------------------------------------------
CS    DATA PLG/-2.4562334077563243311E+01,2.3642701335621505212E+02,
CS   1         -5.4989956895857911039E+02,3.5687548468071500413E+02/
CS    DATA QLG/-3.5553900764052419184E+01,1.9400230218539473193E+02,
CS   1         -3.3442903192607538956E+02,1.7843774234035750207E+02/
      DATA PLG/-2.4562334077563243311D+01,2.3642701335621505212D+02,
     1         -5.4989956895857911039D+02,3.5687548468071500413D+02/
      DATA QLG/-3.5553900764052419184D+01,1.9400230218539473193D+02,
     1         -3.3442903192607538956D+02,1.7843774234035750207D+02/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C  J0(X) / (X**2 - XJ0**2),  XSMALL  <  |X|  <=  4.0
C--------------------------------------------------------------------
CS    DATA PJ0/6.6302997904833794242E+06,-6.2140700423540120665E+08,
CS   1         2.7282507878605942706E+10,-4.1298668500990866786E+11,
CS   2        -1.2117036164593528341E-01, 1.0344222815443188943E+02,
CS   3        -3.6629814655107086448E+04/
CS    DATA QJ0/4.5612696224219938200E+05, 1.3985097372263433271E+08,
CS   1         2.6328198300859648632E+10, 2.3883787996332290397E+12,
CS   2         9.3614022392337710626E+02/
      DATA PJ0/6.6302997904833794242D+06,-6.2140700423540120665D+08,
     1         2.7282507878605942706D+10,-4.1298668500990866786D+11,
     2        -1.2117036164593528341D-01, 1.0344222815443188943D+02,
     3        -3.6629814655107086448D+04/
      DATA QJ0/4.5612696224219938200D+05, 1.3985097372263433271D+08,
     1         2.6328198300859648632D+10, 2.3883787996332290397D+12,
     2         9.3614022392337710626D+02/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C  J0(X) / (X**2 - XJ1**2),  4.0  <  |X|  <=  8.0
C-------------------------------------------------------------------
CS    DATA PJ1/4.4176707025325087628E+03, 1.1725046279757103576E+04,
CS   1         1.0341910641583726701E+04,-7.2879702464464618998E+03,
CS   2        -1.2254078161378989535E+04,-1.8319397969392084011E+03,
CS   3         4.8591703355916499363E+01, 7.4321196680624245801E+02/
CS    DATA QJ1/3.3307310774649071172E+02,-2.9458766545509337327E+03,
CS   1         1.8680990008359188352E+04,-8.4055062591169562211E+04,
CS   2         2.4599102262586308984E+05,-3.5783478026152301072E+05,
CS   3        -2.5258076240801555057E+01/
      DATA PJ1/4.4176707025325087628D+03, 1.1725046279757103576D+04,
     1         1.0341910641583726701D+04,-7.2879702464464618998D+03,
     2        -1.2254078161378989535D+04,-1.8319397969392084011D+03,
     3         4.8591703355916499363D+01, 7.4321196680624245801D+02/
      DATA QJ1/3.3307310774649071172D+02,-2.9458766545509337327D+03,
     1         1.8680990008359188352D+04,-8.4055062591169562211D+04,
     2         2.4599102262586308984D+05,-3.5783478026152301072D+05,
     3        -2.5258076240801555057D+01/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C    (Y0(X) - 2 LN(X/XY0) J0(X)) / (X**2 - XY0**2),
C        XSMALL  <  |X|  <=  3.0
C--------------------------------------------------------------------
CS    DATA PY0/1.0102532948020907590E+04,-2.1287548474401797963E+06,
CS   1         2.0422274357376619816E+08,-8.3716255451260504098E+09,
CS   2         1.0723538782003176831E+11,-1.8402381979244993524E+01/
CS    DATA QY0/6.6475986689240190091E+02, 2.3889393209447253406E+05,
CS   1         5.5662956624278251596E+07, 8.1617187777290363573E+09,
CS   2         5.8873865738997033405E+11/
      DATA PY0/1.0102532948020907590D+04,-2.1287548474401797963D+06,
     1         2.0422274357376619816D+08,-8.3716255451260504098D+09,
     2         1.0723538782003176831D+11,-1.8402381979244993524D+01/
      DATA QY0/6.6475986689240190091D+02, 2.3889393209447253406D+05,
     1         5.5662956624278251596D+07, 8.1617187777290363573D+09,
     2         5.8873865738997033405D+11/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C    (Y0(X) - 2 LN(X/XY1) J0(X)) / (X**2 - XY1**2),
C        3.0  <  |X|  <=  5.5
C--------------------------------------------------------------------
CS    DATA PY1/-1.4566865832663635920E+04, 4.6905288611678631510E+06,
CS   1         -6.9590439394619619534E+08, 4.3600098638603061642E+10,
CS   2         -5.5107435206722644429E+11,-2.2213976967566192242E+13,
CS   3          1.7427031242901594547E+01/
CS    DATA QY1/ 8.3030857612070288823E+02, 4.0669982352539552018E+05,
CS   1          1.3960202770986831075E+08, 3.4015103849971240096E+10,
CS   2          5.4266824419412347550E+12, 4.3386146580707264428E+14/
      DATA PY1/-1.4566865832663635920D+04, 4.6905288611678631510D+06,
     1         -6.9590439394619619534D+08, 4.3600098638603061642D+10,
     2         -5.5107435206722644429D+11,-2.2213976967566192242D+13,
     3          1.7427031242901594547D+01/
      DATA QY1/ 8.3030857612070288823D+02, 4.0669982352539552018D+05,
     1          1.3960202770986831075D+08, 3.4015103849971240096D+10,
     2          5.4266824419412347550D+12, 4.3386146580707264428D+14/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C    (Y0(X) - 2 LN(X/XY2) J0(X)) / (X**2 - XY2**2),
C        5.5  <  |X|  <=  8.0
C--------------------------------------------------------------------
CS    DATA PY2/ 2.1363534169313901632E+04,-1.0085539923498211426E+07,
CS   1          2.1958827170518100757E+09,-1.9363051266772083678E+11,
CS   2         -1.2829912364088687306E+11, 6.7016641869173237784E+14,
CS   3         -8.0728726905150210443E+15,-1.7439661319197499338E+01/
CS    DATA QY2/ 8.7903362168128450017E+02, 5.3924739209768057030E+05,
CS   1          2.4727219475672302327E+08, 8.6926121104209825246E+10,
CS   2          2.2598377924042897629E+13, 3.9272425569640309819E+15,
CS   3          3.4563724628846457519E+17/
      DATA PY2/ 2.1363534169313901632D+04,-1.0085539923498211426D+07,
     1          2.1958827170518100757D+09,-1.9363051266772083678D+11,
     2         -1.2829912364088687306D+11, 6.7016641869173237784D+14,
     3         -8.0728726905150210443D+15,-1.7439661319197499338D+01/
      DATA QY2/ 8.7903362168128450017D+02, 5.3924739209768057030D+05,
     1          2.4727219475672302327D+08, 8.6926121104209825246D+10,
     2          2.2598377924042897629D+13, 3.9272425569640309819D+15,
     3          3.4563724628846457519D+17/
C-------------------------------------------------------------------
C  Coefficients for Hart,s approximation,  |X| > 8.0
C-------------------------------------------------------------------
CS    DATA P0/3.4806486443249270347E+03, 2.1170523380864944322E+04,
CS   1        4.1345386639580765797E+04, 2.2779090197304684302E+04,
CS   2        8.8961548424210455236E-01, 1.5376201909008354296E+02/
CS    DATA Q0/3.5028735138235608207E+03, 2.1215350561880115730E+04,
CS   1        4.1370412495510416640E+04, 2.2779090197304684318E+04,
CS   2        1.5711159858080893649E+02/
CS    DATA P1/-2.2300261666214198472E+01,-1.1183429920482737611E+02,
CS   1        -1.8591953644342993800E+02,-8.9226600200800094098E+01,
CS   2        -8.8033303048680751817E-03,-1.2441026745835638459E+00/
CS    DATA Q1/1.4887231232283756582E+03, 7.2642780169211018836E+03,
CS   1        1.1951131543434613647E+04, 5.7105024128512061905E+03,
CS   2        9.0593769594993125859E+01/
      DATA P0/3.4806486443249270347D+03, 2.1170523380864944322D+04,
     1        4.1345386639580765797D+04, 2.2779090197304684302D+04,
     2        8.8961548424210455236D-01, 1.5376201909008354296D+02/
      DATA Q0/3.5028735138235608207D+03, 2.1215350561880115730D+04,
     1        4.1370412495510416640D+04, 2.2779090197304684318D+04,
     2        1.5711159858080893649D+02/
      DATA P1/-2.2300261666214198472D+01,-1.1183429920482737611D+02,
     1        -1.8591953644342993800D+02,-8.9226600200800094098D+01,
     2        -8.8033303048680751817D-03,-1.2441026745835638459D+00/
      DATA Q1/1.4887231232283756582D+03, 7.2642780169211018836D+03,
     1        1.1951131543434613647D+04, 5.7105024128512061905D+03,
     2        9.0593769594993125859D+01/
C-------------------------------------------------------------------
C  Check for error conditions
C-------------------------------------------------------------------
      AX = ABS(ARG)
      IF ((JINT .EQ. 1) .AND. (ARG .LE. ZERO)) THEN
            RESULT = -XINF
            GO TO 2000
         ELSE IF (AX .GT. XMAX) THEN
            RESULT = ZERO
            GO TO 2000
      END IF
      IF (AX .GT. EIGHT) GO TO 800
      IF (AX .LE. XSMALL) THEN
         IF (JINT .EQ. 0) THEN
               RESULT = ONE
            ELSE
               RESULT = PI2 * (LOG(AX) + CONS)
         END IF
         GO TO 2000
      END IF
C-------------------------------------------------------------------
C  Calculate J0 for appropriate interval, preserving
C     accuracy near the zero of J0
C-------------------------------------------------------------------
      ZSQ = AX * AX
      IF (AX .LE. FOUR) THEN
            XNUM = (PJ0(5) * ZSQ + PJ0(6)) * ZSQ + PJ0(7)
            XDEN = ZSQ + QJ0(5)
            DO 50 I = 1, 4
               XNUM = XNUM * ZSQ + PJ0(I)
               XDEN = XDEN * ZSQ + QJ0(I)
   50       CONTINUE
            PROD = ((AX - XJ01/TWO56) - XJ02) * (AX + XJ0)
         ELSE
            WSQ = ONE - ZSQ / SIXTY4
            XNUM = PJ1(7) * WSQ + PJ1(8)
            XDEN = WSQ + QJ1(7)
            DO 220 I = 1, 6
               XNUM = XNUM * WSQ + PJ1(I)
               XDEN = XDEN * WSQ + QJ1(I)
  220       CONTINUE
            PROD = (AX + XJ1) * ((AX - XJ11/TWO56) - XJ12)
      END IF
      RESULT = PROD * XNUM / XDEN
      IF (JINT .EQ. 0) GO TO 2000
C-------------------------------------------------------------------
C  Calculate Y0.  First find  RESJ = pi/2 ln(x/xn) J0(x),
C    where xn is a zero of Y0
C-------------------------------------------------------------------
      IF (AX .LE. THREE) THEN
            UP = (AX-XY01/TWO56)-XY02
            XY = XY0
         ELSE IF (AX .LE. FIVE5) THEN
            UP = (AX-XY11/TWO56)-XY12
            XY = XY1
         ELSE
            UP = (AX-XY21/TWO56)-XY22
            XY = XY2
      END IF
      DOWN = AX + XY
      IF (ABS(UP) .LT. P17*DOWN) THEN
            W = UP/DOWN
            WSQ = W*W
            XNUM = PLG(1)
            XDEN = WSQ + QLG(1)
            DO 320 I = 2, 4
               XNUM = XNUM*WSQ + PLG(I)
               XDEN = XDEN*WSQ + QLG(I)
  320       CONTINUE
            RESJ = PI2 * RESULT * W * XNUM/XDEN
         ELSE
            RESJ = PI2 * RESULT * LOG(AX/XY)
      END IF
C-------------------------------------------------------------------
C  Now calculate Y0 for appropriate interval, preserving
C     accuracy near the zero of Y0
C-------------------------------------------------------------------
      IF (AX .LE. THREE) THEN
            XNUM = PY0(6) * ZSQ + PY0(1)
            XDEN = ZSQ + QY0(1)
            DO 340 I = 2, 5
               XNUM = XNUM * ZSQ + PY0(I)
               XDEN = XDEN * ZSQ + QY0(I)
  340       CONTINUE
         ELSE IF (AX .LE. FIVE5) THEN
            XNUM = PY1(7) * ZSQ + PY1(1)
            XDEN = ZSQ + QY1(1)
            DO 360 I = 2, 6
               XNUM = XNUM * ZSQ + PY1(I)
               XDEN = XDEN * ZSQ + QY1(I)
  360       CONTINUE
         ELSE
            XNUM = PY2(8) * ZSQ + PY2(1)
            XDEN = ZSQ + QY2(1)
            DO 380 I = 2, 7
               XNUM = XNUM * ZSQ + PY2(I)
               XDEN = XDEN * ZSQ + QY2(I)
  380       CONTINUE
      END IF
      RESULT = RESJ + UP * DOWN * XNUM / XDEN
      GO TO 2000
C-------------------------------------------------------------------
C  Calculate J0 or Y0 for |ARG|  >  8.0
C-------------------------------------------------------------------
  800 Z = EIGHT / AX
      W = AX / TWOPI
      W = AINT(W) + ONEOV8
      W = (AX - W * TWOPI1) - W * TWOPI2
      ZSQ = Z * Z
      XNUM = P0(5) * ZSQ + P0(6)
      XDEN = ZSQ + Q0(5)
      UP = P1(5) * ZSQ + P1(6)
      DOWN = ZSQ + Q1(5)
      DO 850 I = 1, 4
         XNUM = XNUM * ZSQ + P0(I)
         XDEN = XDEN * ZSQ + Q0(I)
         UP = UP * ZSQ + P1(I)
         DOWN = DOWN * ZSQ + Q1(I)
  850 CONTINUE
      R0 = XNUM / XDEN
      R1 = UP / DOWN
      IF (JINT .EQ. 0) THEN
            RESULT = SQRT(PI2/AX) * (R0*COS(W) - Z*R1*SIN(W))
         ELSE
            RESULT = SQRT(PI2/AX) * (R0*SIN(W) + Z*R1*COS(W))
      END IF
 2000 RETURN
C---------- Last line of CALJY0 ----------
      END
*DECK J4SAVE
      FUNCTION J4SAVE (IWHICH, IVALUE, ISET)
C***BEGIN PROLOGUE  J4SAVE
C***SUBSIDIARY
C***PURPOSE  Save or recall global variables needed by error
C            handling routines.
C***LIBRARY   SLATEC (XERROR)
C***TYPE      INTEGER (J4SAVE-I)
C***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        J4SAVE saves and recalls several global variables needed
C        by the library error handling routines.
C
C     Description of Parameters
C      --Input--
C        IWHICH - Index of item desired.
C                = 1 Refers to current error number.
C                = 2 Refers to current error control flag.
C                = 3 Refers to current unit number to which error
C                    messages are to be sent.  (0 means use standard.)
C                = 4 Refers to the maximum number of times any
C                     message is to be printed (as set by XERMAX).
C                = 5 Refers to the total number of units to which
C                     each error message is to be written.
C                = 6 Refers to the 2nd unit for error messages
C                = 7 Refers to the 3rd unit for error messages
C                = 8 Refers to the 4th unit for error messages
C                = 9 Refers to the 5th unit for error messages
C        IVALUE - The value to be set for the IWHICH-th parameter,
C                 if ISET is .TRUE. .
C        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
C                 given the value, IVALUE.  If ISET=.FALSE., the
C                 IWHICH-th parameter will be unchanged, and IVALUE
C                 is a dummy parameter.
C      --Output--
C        The (old) value of the IWHICH-th parameter will be returned
C        in the function value, J4SAVE.
C
C***SEE ALSO  XERMSG
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900205  Minor modifications to prologue.  (WRB)
C   900402  Added TYPE section.  (WRB)
C   910411  Added KEYWORDS section.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  J4SAVE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
C***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
      SUBROUTINE MAGMIE(NANGL,CXEPS,CXMU,A,LAMBDA,
     &                  QABS,QEXT,QSCA,G,THETA,S11,S12,S33,S34,POL)
      IMPLICIT NONE 

c Arguments:

      INTEGER NANGL
      DOUBLE PRECISION A,G,LAMBDA,QABS,QEXT,QSCA
      DOUBLE PRECISION
     &   POL(1:NANGL),
     &   S11(1:NANGL),
     &   S12(1:NANGL),
     &   S33(1:NANGL),
     &   S34(1:NANGL),
     &   THETA(1:NANGL)
      DOUBLE COMPLEX CXEPS,CXMU

c local variables:

      INTEGER NANGLMX
      PARAMETER(NANGLMX=255)
      DOUBLE PRECISION ONE,ZERO
      PARAMETER(ONE=1.D0,ZERO=0.D0) 

      DOUBLE COMPLEX S1(NANGLMX),S2(NANGLMX),SX1(NANGLMX),SX2(NANGLMX)
      DOUBLE COMPLEX CXEPSN,CXMUN

      DOUBLE PRECISION I1,I2,INTEN,XX,ANGINC
      DOUBLE PRECISION ISX1,ISX2,INTSX,DGPLSX,ANGLE,DEGPOL
      DOUBLE PRECISION PI
      DOUBLE PRECISION S11NOR
      DOUBLE PRECISION SX11(NANGLMX),SX12(NANGLMX),SX33(NANGLMX)
      DOUBLE PRECISION SX34(NANGLMX),SXPOL(NANGLMX)
      DOUBLE PRECISION QBAC
      DOUBLE PRECISION QESX,QSSX,QASX,QBSX

      INTEGER SIGN,MIOR,I,NX,NMANG1,NMANG2,NIOR,NXX 

      LOGICAL XFLAG

c *********************************************************************
c Subroutine MAGMIE

c given:
c     NANGL = number of angles from 0 to 180 (including 0 and 180)
c             (should be ODD) so that 90 is included.

c based on program MSMAIN written by M.E. Milham
c and coded and tested by Thomas Wriedt (see below)
c
c     Program coded and tested by 
c     Thomas Wriedt, thw@iwt.uni-bremen.de, www.ScattPort.org
c     Bremen, Germany, 9. March 2010
c     
c     according to the report
c     Merril E. Milham, Electromagnetic scattering by magnetic spheres: 
c     theory and  algorithms, 
c     Edgewood Research Development and Engineering Center, 
c     ERDEC-TR-207,  Oct. 1994, ADA289798.
c
c     The program is available via 
c     www.scattport.org/index.php/programs-menu/
c     mie-type-codes-menu/70-mie-code-magnetic-sphere
c
c     Compiler: Intel(R) Fortran Compiler
c 
c     Program changes: 
c     - missing SLATEC subroutines added
c     - zabs and zsqrt replayed by abs and sqrt 
c     - parameter input via input.dat file
c     - pages numbers refer to the report
c     - changes are marked by cthw
c
c *************************************************************************
c
c Program msmain computes the electromagnetic scattering for magnetic
c spheres. If the sphere is small a comparison is made between the
c results of a rigorous, complete calculation and an approximate
c calculation valid for small spheres.
c
c     Merrill Milham    >>> version 1.0 <<<     JANUARY 1994
c
c inputs:
c               ALL INPUT IS IN LIST DIRECTED FORMAT
c
c     line #1
c               flag = 'r' for refractive index data or
c                      'p' for permittivity data (character*l)
c               mior = number of complex indexes or
c                      complex permittivities to be read (integer)
c                 mp = complex index or permittivity values
c                       to be read, up to 10 values (complex*16)
c     line #2
c                muu = complex permeability values to be read
c                      mior values are required (complex*16)
c     line #3
c                 nx = number of size parameters to be read (integer)
c                  x = size parameter values to be read (real*8)
c     line #4 anginc = angular increment in degrees which is
c                      added to zero to produce the angles at
c                      which the scattering is to be produced. 
c 		       Values of anginc are restricted to those
c                      for which mod(90,anginc).eq.zero. (real*8)
C
c output:
c           For each angle the following quantities are given:
c                   The complex scattered field amplitudes, the
c                   scattered intensity, the degree of polarization,
c                   the Mueller matrix elements, and the polarization.
C
c           Efficiency factors for extinction, scattering, absorption, and
c           backscatter, and the asymmetry factor.
C
c subroutines used:
C
c            magsph - returns computed scattering quantities for a
c                     homgeneous magnetic sphere or homogeneous
c                     nonmagnetic sphere if the permeability = (1,O)
C
c            amuelr - returns elements of the Mueller matrix for a
c                     spherical scatterer
C
c            msphsx - returns approximate scattering quantities for a
c                     small homogeneous, magnetic sphere or a small,
c                     homogeneous, nonmagnetic sphere if the
c                     permeability = (1,0)
c history
c 2011.03.07 (BTD) downloaded from scattport library maintained by Thomas Wriedt
c              http://www.scattport.org/index.php/programs-menu/
c                     mie-type-codes-menu/70-mie-code-magnetic-sphere
c 2011.03.08 (BTD) cosmetic changes to code to make it more readable
c 2011.03.12 (BTD) 
c end history             
c **************************************************************************
c-----------------------------------------------------------------------
c*** diagnostic
c      write(0,*)'magmie ckpt 1'
c      write(0,*)' nangl=',nangl
c      write(0,*)' cxeps=',cxeps
c      write(0,*)' cxmu=',cxmu
c      write(0,*)' a=',a
c      write(0,*)' lambda=',lambda
c***

      PI=4.D0*ATAN(1.D0)

c sanity checks -----------------------------------------------------

      IF(NANGL.GT.NANGLMX)THEN
         WRITE(0,*)'magmie fatal error: NANGL=',NANGL,' > NANGLMX=',
     &             NANGLMX
         WRITE(0,*)' need to recompile with larger NANGLMX'
         STOP
      ENDIF
      IF(MOD(NANGL,2).NE.1)THEN
         WRITE(0,*)'magmie fatal error: NANGL=',NANGL,
     &             ' must be ODD'
         STOP
      ENDIF

c end sanity checks -------------------------------------------------

      NMANG1=NANGL
      NMANG2=(NANGL+1)/2

      DO I=1,NMANG1
         THETA(I)=DBLE((I-1)*180)/DBLE(NANGL-1)
      ENDDO

      SIGN=1

      XX=2.*PI*A/LAMBDA

      XFLAG=XX*ABS(SQRT(CXMU*CXEPS)).LT.ONE

c convention here is that Im(epsilon) .le. 0
c                         Im(mu)      .le. 0

      CXEPSN=CXEPS
      IF(IMAG(CXEPS).GT.0.)CXEPSN=CONJG(CXEPS)
      CXMUN=CXMU
      IF(IMAG(CXMU).GT.0.)CXMUN=CONJG(CXMU)

c*** diagnostic
c      write(0,*)'magmie ckpt 100'
c      write(0,*)' xx=',xx
c      write(0,*)' cxepsn=',cxepsn
c      write(0,*)' cxmun=',cxmun
c      write(0,*)' nmang2=',nmang2
c      write(0,*)' theta=',theta
c***
      CALL MAGSPH(XX,CXEPSN,CXMUN,NMANG2,THETA,QEXT,QSCA,QBAC,G,S1,S2)
c*** diagnostic
c      write(0,*)'magmie ckpt 101 (returned from magsph)'
c      write(0,*)' QEXT=',QEXT
c***
      QABS=QEXT-QSCA

      S11NOR=0.D0
      CALL AMUELR(S1,S2,THETA,NMANG1,SIGN,S11NOR,S11,S12,S33,S34,POL)

      IF(XFLAG)THEN
c*** diagnostic
c         write(0,*)'magmie ckpt 110: call msphsx'
c***
         CALL MSPHSX(XX,CXMUN,CXEPSN,NMANG2,THETA,QESX,QSSX,QASX,QBSX,
     &               SX1,SX2)
c*** diagnostic
c         write(0,*)'magmie ckpt 111'
c         write(0,*)' qesx=',qesx
c***
         S11NOR=0.D0
         CALL AMUELR(SX1,SX2,THETA,NMANG1,SIGN,S11NOR,SX11,SX12,
     &               SX33,SX34,SXPOL)
      ELSE
         CONTINUE
      ENDIF

      DO I=1,NMANG1
         ANGLE=THETA(I)

         I1=(DBLE(S1(I)))**2+(DIMAG(S1(I)))**2
         I2=(DBLE(S2(I)))**2+(DIMAG(S2(I)))**2
         INTEN=0.5*(I1+I2)
         IF(.NOT.(I1.EQ.ZERO.AND.I2.EQ.ZERO))THEN
            DEGPOL=(I2-I1)/(I2+I1)
         ELSE
            WRITE(*,*)
            WRITE(*,*)'degree of polarization undefined for ',
     &                ANGLE,' deg'
            WRITE(*,*)
            DEGPOL=ZERO
         ENDIF

         IF(XFLAG)THEN
            ISX1=(DBLE(SX1(I)))**2+(DIMAG(SX1(I)))**2
            ISX2=(DBLE(SX2(I)))**2+(DIMAG(SX2(I)))**2
            INTSX=0.5*(ISX1+ISX2)
            IF(.NOT.(ISX1.EQ.ZERO.AND.ISX2.EQ.ZERO))THEN
               DGPLSX=(ISX2-ISX1)/(ISX2+ISX1)
            ELSE
               WRITE(*,*)
               WRITE(*,*)'degree of polarization undefined for ',
     &                   ANGLE,' deg.sx'
               WRITE(*,*)
               DGPLSX=ZERO
            ENDIF
         ELSE
            CONTINUE
         ENDIF 
      ENDDO
c*** diagnostic
c      write(0,*)'magmie ckpt 900'
c      write(0,*)' qext=',qext
c***
      RETURN
      END
      SUBROUTINE MAGSPH(X,EPS,MU,NUMANG,THETA,
     &                  QEXT,QSCA,QBAC,G,S1,S2)
c
c ****************************************************
c
c Subroutine magsph computes the scattering cross sections and angular
c scattering from a magnetic sphere. If the number of scattering angles
c is set to zero, only the cross sections (efficiencies) are returned.
c
c     Merrill Milham        >>> version 2.0 <<<            SEPT 1993
c
c     Inputs :
c               x = size parameter of the sphere           (real*8)
c             eps = complex permittivity: epsr -i*epsi     (complex*16)
c              mu = complex permeability: mur - i*mui      (complex*16)
c          numang = number of scattering angles            (integer)
c                   between 0 & 90 deg.
c           theta = scattering angles in degrees           (real*8)
c                theta(i) are entered between 0 & 90 deg.
c                theta must increase monotonically. Results for
c                supplementary angles (180 deg. - theta(i)) are
c                also returned.
C
c     Outputs :
c           qext = extinction efficiency
c           qsca = scattering efficiency
c           qbac = backscatter efficiency
c               g = asymmetry factor
c              s1 = scattered amplitude
c              s2 = scattered amplitude
C
c Subroutines used:
c
c           zbjy - returns one-half integer order J & Y Bessel functions
c
c References:
c M. Kerker, D.4. Wang, and C. L. Giles, "Electromagnetic scattering
c by magnetic spheres," J. Opt. Soc. Am., 73, 765-767 (1983).
c
c D. E. Amos, "Algorithm 644:A portable package for Bessel functions of
c a complex argument and nonnegative order," ACM Trans. on Math.
c Software, 12, 265-273 (1986).
c
c M. Abramowitz and I. A. Stegun, "Handbook of Mathematical Functions,"
c NBS Applied Math. Series 55, US Dept. of Commerce, Washington, DC
c (1955).
c
c W. J. Wiscombe,"Mie Scattering Calculations: Advances in Technique
c
c page 25 
c
c and Fast, Vector-Speed Computer Codes," NCAR Tech. Note, NCAR/TN-14O+STR
c (1979)
c
c .........................................................................
c
      IMPLICIT NONE
c
      REAL*8 X
      COMPLEX*16 EPS,MU
      INTEGER NUMANG
      REAL*8 THETA(1) 
*
      REAL*8 QEXT,QSCA,QBAC,G
      COMPLEX*16 S1(1),S2(1)
c
      INTEGER AL,NANGL,NANGL2
      REAL*8 THIRD
      PARAMETER (THIRD=1.D0/3.D0,AL=5100,NANGL=255,NANGL2=(NANGL+1)/2)
      COMPLEX*16 SP(NANGL2),SM(NANGL2),SPS(NANGL2),SMS(NANGL2) 
*
      COMPLEX*16 M,MC1,MXI,S,T,U,V,AN,BN,XP
      REAL*8 XI,DN,DNN,RN,TNP1,THETAN
      REAL*8 XMU(NANGL),PI(NANGL),PIL(NANGL),TAU(NANGL) 
*
      REAL*8 BJR(AL),BYR(AL)
      REAL*8 CJR(AL),CJI(AL),CYR(AL),CYI(AL)
      REAL*8 BJN,BYN,BJL
      REAL*8 SC,CA,T1,T2,T3,T4
      COMPLEX*16 CJN,CYN,CJL,CYL,H2N,H2L
      COMPLEX*16 ANL,BNL,BS,ANP,BNP,ABP,ABM,ANPM,BNPM
      COMPLEX*16 ZT1,ZT2 
*
      INTEGER KSTOP,K,MM,N,J,J2,J3 
*
      REAL*8 CPI,ZERO,ONE,TWO,RAD,HALF,FNU
      COMPLEX*16 CDBLEI,CDBLE1,CDBLE0
      PARAMETER (CPI=3.1415926535897932384D0,ZERO=0.D0,ONE=1.D0)
      PARAMETER (TWO=2.D0,RAD=CPI/180.D0,HALF=0.5D0,FNU=HALF)
      PARAMETER (CDBLEI=(0.D0,1.D0),CDBLE1=(1.D0,0.D0))
      PARAMETER (CDBLE0=(0.D0,0.D0))
c
      KSTOP=IDINT(X+4.D0*X**THIRD+4.D0) 
*
      IF(KSTOP.LE.AL) THEN
      CONTINUE
                      ELSE
      PRINT*,'magsph arrays too small: kstop =',KSTOP,' al =',AL
      STOP
      END IF
c
c page 26 
c
      IF(NUMANG.EQ.0) THEN
      S1(1)=CDBLE0
      S2(1)=CDBLE0
                      ELSE 
*
            IF(NUMANG.LE.NANGL2) THEN
                CONTINUE
                                 ELSE
        PRINT*,NUMANG,'scattering angles input: only',NANGL2,' allowed'
            STOP
           ENDIF 
*
      DO 100 N=1,NUMANG
        
      THETAN=DABS(THETA(N))
      THETA(N)=THETAN
*
           IF(THETAN.LE.90.D0) THEN
               CONTINUE
                                   ELSE
      PRINT*,'theta(',n,')=',THETAN,'scattering angles must be < 90 deg'
            STOP
            END IF 
*
      THETAN=RAD*THETAN
      XMU(N)=DCOS(THETAN) 
*
      SP(N)=CDBLE0
      SM(N)=CDBLE0
      SPS(N)=CDBLE0
      SMS(N)=CDBLE0
      PI(N)=HALF
      PIL(N)=ZERO 
*
100   CONTINUE 
*
      END IF
c
cthw      m=zsqrt(mu*eps)
      M=SQRT(MU*EPS)

      MC1=M/MU 
*
      XI=ONE/X
      MXI=XI/M
c
      CALL ZBJY(X,M,KSTOP,FNU,BJR,BYR,CJR,CJI)
C
      BJL=BJR(1)
      CJL=DCMPLX(CJR(1),CJI(1))
cthw      cyl=dcmplx(cyr(1),cyi(1))
      H2L=DCMPLX(BJR(1),-BYR(1))
*
      QEXT=ZERO
      QSCA=ZERO
c
c page 27 
c
      BS=CDBLE0
      G=ZERO
*
      ANL=CDBLE0
      BNL=CDBLE0
      DN=ONE
      RN=ONE
cthw  tnp1-one 
      TNP1=ONE
      MM=1
c
      DO 300 K=2,KSTOP
*
      TNP1=TNP1+TWO
      T1=DN-RN
      CA=ONE+RN
      SC=RN
*
      DNN=DN+ONE
      RN=ONE/DNN
*
      SC=SC+RN
* 
      BJN=BJR(K)
      BYN=BYR(K)
      CJN=DCMPLX(CJR(K),CJI(K))
cthw      cyn=dcmplx(cyr(k),cyi(k))
      H2N=DCMPLX(BJN,-BYN)
*      
      XP=DN*MXI
      S=CJL-XP*CJN
      U=S*H2N
      S=S*BJN
*
      XP=DN*DCMPLX(XI,ZERO)
      T=CJN*(BJL-XP*BJN)
      V=CJN*(H2L-XP*H2N)
*      
      AN=(S-MC1*T)/(U-MC1*V)
      BN=(MC1*S-T)/(MC1*U-V)
      ABP=AN+BN
      ABM=AN-BN
*
      ZT1=DCONJG(AN)
      ZT2=DCONJG(BN)
      QEXT=QEXT+TNP1*DBLE(ABP)
      QSCA=QSCA+TNP1*(AN*ZT1+BN*ZT2)
      BS=BS-(DN+HALF)*MM*ABM
      G=G+T1*DBLE(ANL*ZT1+BNL*ZT2)+SC*DBLE(AN*ZT2)
*
      IF(NUMANG.EQ.0) THEN
        CONTINUE
                      ELSE
c
c page 28 
c
      ANP=SC*ABP
      BNP=SC*ABM
      ANPM=MM*ANP
      BNPM=MM*BNP
*
      DO 375 J=1,NUMANG
      T1=XMU(J)*PI(J)
      T4=T1-PIL(J)
      TAU(J)=DN*T4-PIL(J)
      T2=PI(J)+TAU(J)
      T3=PI(J)-TAU(J) 
*
      SP(J)=SP(J)+ANP*T2
      SMS(J)=SMS(J)+BNPM*T2
      SM(J)=SM(J)+BNP*T3
      SPS(J)=SPS(J)+ANPM*T3 
*
      PIL(J)=PI(J)
      PI(J)=T1+CA*T4 
*
375   CONTINUE
      END IF 
*
      DN=DNN
      MM=-MM
      ANL=AN
      BNL=BN 
*
      BJL=BJN
      CJL=CJN
cthw      cyl=cyn
      H2L=H2N 
*
300   CONTINUE
C
      IF(NUMANG.EQ.0) THEN
      CONTINUE
                        ELSE
      	J2=2*NUMANG
      	DO 500 J=1,NUMANG
      	J3=J2-J
      	S1(J)=SP(J)+SM(J)
      	S2(J)=SP(J)-SM(J)
      	S1(J3)=SPS(J)+SMS(J)
      	S2(J3)=SPS(J)-SMS(J)
500   	CONTINUE
*
      END IF
C
      XI=TWO*XI*XI
      QEXT=XI*QEXT
      QSCA=XI*QSCA
c
c page 29
c
      XI=TWO*XI
      QBAC=XI*BS*DCONJG(BS)
      G=XI/QSCA*G
c
      RETURN
c
      END
      SUBROUTINE MIEV0( XX, CREFIN, PERFCT, MIMCUT, ANYANG, NUMANG, XMU,
     &                  NMOM, IPOLZN, MOMDIM, PRNT, QEXT, QSCA, GQSC,
     &                  PMOM, SFORW, SBACK, S1, S2, TFORW, TBACK,
     &                  SPIKE )

c Subroutine MIEV0: double precision version

c Given:

c   XX     = Mie size parameter = 2*pi*radius/lambda,
c            where lambda = wavelength in surrounding medium
c   CREFIN = complex refractive index
c            (imaginary part can be given as + or -, but calculation
c             will be carried out for an absorbing material)
c   PERFCT = .FALSE. to calculate S1,S2 for finite refractive index
c            .TRUE. to calculate S1,S2 for perfect conductor
c                   scattering amplitudes will be returned using
c                   * Bohren & Huffman convention if Im(CREFIN) > 0,
c                   * van de Hulst convention if Im(CREFIN) <= 0
c   MIMCUT = positive value below which imaginary part of refractive
c            index is regarded as zero (computation proceeds faster for
c            zero imaginary index).  
c            0 for standard use.
c   ANYANG = .TRUE. if arbitrary scattering angles are input through XMU
c            .FALSE. if angles are monotone increasing and mirror symmetric
c                    around 90 degrees (computation is especially efficient
c                    for this case)
c   NUMANG = number of angles at which S1 and S2 are to be evaluated
c            0 to skip calculation of S1 and S2
c            NUMANG must not exceed parameter MAXANG !
c   XMU(J) = cosines of scattering angles J=1 - NUMANG
c            at which S1,S2 are to be evaluated
c            If ANYANG=.FALSE., these angles must be monotone increasing
c                               and mirror-symmetric about 90 deg.
c                               if ANYANG=.FALSE. and NUMANG is odd,
c                               90 deg. must be among the angles
c   NMOM   = highest Legendre moment PMOM to calculate, numbering from zero
c          = 0 to suppress computation of Legendre moments
c   IPOLZN > 0 : if NMON > 0, compute Legendre moments PMOM for Mueller 
c                matrix elements determined by digits of IPOLZN:
c                   1 refers to M1
c                   2           M2
c                   3           S21
c                   4           D21
c                e.g., if IPOLZN=14, PMOM will be evaluated only for M1
c                and D21
c          = 0 : if NMOM > 0, compute Legendre moments PMOM only for the 
c                unpolarized unnormalized phase function (S11 ?)
c          < 0 : if NMOM > 0, compute Legendre moments PMOM for the
c                Sekera phase quantities determined by the digits of
c                ABS(IPOLZN):
c                    1 refers to R1
c                    2           R2
c                    3           R3
c                    4           R4
c                e.g., if IPOLZN=-14, PMOM will be evaluated only for
c                R1 and R4
c   MOMDIM = DIMENSIONing information for PMOM(0:MOMDIM,*)
c            Note that PMOM must be DIMENSIONed in the calling program,
c            with second dimension the larger of unity and the highest
c            digit in IPOLZN
c   PRNT(1)= .TRUE. to print S1,S2,|S1|^2,|S2|^2, and polarization
c            .FALSE. to suppress printing S1,S2,|S1|^2,|S2|^2,pol
c   PRNT(2)= .TRUE. to print all output variables other than
c                            S1,S2,|S1|^2,|S2|^2,pol
c            .FALSE. to suppress printing these additional variables

c Returns:
c   QEXT   = C_ext/(pi*a^2)
c   QSCA   = C_sca/(pi*a^2)
c   GQSC   = <cos(theta)>*QSCA
c   S1(J)  = complex scattering amplitude S_1 for directions J=1-NUMANG
c   S2(J)  = complex scattering amplitude S_2 for directions J=1-NUMANG
c   SFORW  = forward scattering amplitude S_1(theta=0)=S_2(theta=0)
c   SBACK  = complex scattering amplitude S_1(theta=180)=-S_2(theta=180)
c   TFORW(1)=lim [S_2(theta)-cos(theta)*S_1(theta)]/[1-cos^2(theta)]
c            theta->0
c   TFORW(2)=lim [S_1(theta)-cos(theta)*S_2(theta)]/[1-cos^2(theta)]
c            theta->0
c            TFORW(1) and TFORW(2) are required for polarized radiative
c            transfer calculations
c   TBACK(1)=lim [S_2(theta)-cos(theta)*S_1(theta)]/[1-cos^2(theta)]
c            theta->180
c   TBACK(2)=lim [S_1(theta)-cos(theta)*S_2(theta)]/[1-cos^2(theta)]
c            theta->180
c   SPIKE  = (real) magnitude of smallest denominator of either a_n
c            or b_n taken over all terms in the Mie series past
c            N=size paramter XX.  Values of SPIKE < 0.3 signify a ripple
c            spike, since these spikes are produced by abnormally small
c            denominators in the Mie coefficients (normal denominators
c            are of order unity or higher).  Defaults to 1.0 when not on
c            a spike.  Does not identify all resonances.
c   PMOM(N,NP)=real moments N=0-NMOM of unpolarized NP-th phase quantity PQ

c            Moments with M > 2*NTRM are zero, where NTRM = number of
c            terms in Mie series = XX + 4.*XX^(1/3) + 1

c                       inf
c            PQ(mu,NP)= Sum (2m+1) * PMOM(m,NP) * P_m(mu)
c                       m=0

c            where mu = cos(theta)
c                  P_m = m-th Legendre polynomial

c            and the definition of PQ is as follows:
c            IPOLZN > 0:
c                        PQ(mu,1) = |S_1|^2              = M1 = S_11+S_12
c                        PQ(mu,2) = |S_2|^2              = M2 = S_11-S_12
c                        PQ(mu,3) = Re[S_1*conjg(S_2)]   = M3 = S_33
c                        PQ(mu,4) = Im[S_1*conjg(S_2)]   = M4 = S_43
c            IPOLZN = 0:
c                        PQ(mu,1) = (|S_1|^2+|S_2|^2)/2       = S_11
c                                 = "unnormalized phase function"
c            IPOLZN < 0:
c                        PQ(mu,1) = |T_1|^2              = R1
c                        PQ(mu,2) = |T_2|^2              = R2
c                        PQ(mu,3) = Re[T_1*conjg(T_2)]   = R3
c                        PQ(mu,4) = Im[T_1*conjg(T_2)]   = R4
c            note that sign of PQ(mu,4) is as above for convention where
c            Im(refractive index) > 0
c            See Van de Hulst (1957, 1982) for correct formulae for PMOM
c            WARNING: make sure the second dimension of PMOM in the calling
c            program is at least as large as the number of digits in
c            ABS(IPOLZN)
c
c-----------------------------------------------------------------------

c Subroutine MIEV0:
c    Computes Mie scattering and extinction efficiencies; asymmetry
c    factor;  forward- and backscatter amplitude;  scattering
c    amplitudes vs. scattering angle for incident polarization parallel
c    and perpendicular to the plane of scattering;
c    coefficients in the Legendre polynomial expansions of either the
c    unpolarized phase function or the polarized phase matrix;
c    some quantities needed in polarized radiative transfer;  and
c    information about whether or not a resonance has been hit.

c    Input and output variables are described above.
c    Many statements are accompanied by comments referring to 
c    references in MIEV.doc, notably the NCAR Mie report which is now
c    available electronically and which is referred to using the
c    shorthand (Rn), meaning Eq. (n) of the report.

c    CALLING TREE:

c        MIEV0
c            TESTMI
c                TSTBAD
c                MIPRNT
c                ERRMSG
c            CKINMI
c                WRTBAD
c                WRTDIM
c                ERRMSG
c            SMALL1
c            SMALL2
c            ERRMSG
c            BIGA
c                CONFRA
c                    ERRMSG
c            LPCOEF
c                LPCO1T
c                LPCO2T
c                ERRMSG
c            MIPRNT


c      I N T E R N A L   V A R I A B L E S
c      -----------------------------------

c  AN,BN           Mie coefficients a-sub-n, b-sub-n ( Ref. 1, Eq. 16 )

c  ANM1,BNM1       Mie coefficients  a-sub-(n-1),
c                     b-sub-(n-1);  used in GQSC sum

c  ANP             Coeffs. in S+ expansion ( Ref. 2, p. 1507 )
c  BNP             Coeffs. in S- expansion ( Ref. 2, p. 1507 )
c  ANPM            Coeffs. in S+ expansion ( Ref. 2, p. 1507 )
c                     when  MU  is replaced by  - MU
c  BNPM            Coeffs. in S- expansion ( Ref. 2, p. 1507 )
c                     when  MU  is replaced by  - MU

c  CALCMO(K)       TRUE, calculate moments for K-th phase quantity
c                     (derived from IPOLZN)

c  CBIGA(N)        Bessel function ratio A-sub-N (Ref. 2, Eq. 2)
c                     ( COMPLEX version )

c  CDENAN,         (COMPLEX) denominators of An,Bn
c   CDENBN

c  CIOR            Complex index of refraction with negative
c                     imaginary part (Van de Hulst convention)
c  CIORIV          1 / cIoR

c  COEFF           ( 2N + 1 ) / ( N ( N + 1 ) )

c  CSUM1,2         temporary sum variables for TFORW, TBACK

c  FN              Floating point version of loop index for
c                     Mie series summation

c  LITA,LITB(N)    Mie coefficients An, Bn, saved in arrays for
c                     use in calculating Legendre moments PMOM

c  MAXTRM          Max. possible no. of terms in Mie series

c  MM              (-1)^(n+1), where n is Mie series sum index 

c  MIM             Magnitude of imaginary refractive index
c  MRE             Real part of refractive index

c  MAXANG          Max. possible value of input variable NUMANG
c  NANGD2          (NUMANG+1)/2 ( no. of angles in 0-90 deg; ANYANG=F )

c  NOABS           TRUE, sphere non-absorbing (determined by MIMCUT)

c  NP1DN           ( N + 1 ) / N

c  NPQUAN          Highest-numbered phase quantity for which moments are
c                     to be calculated (the largest digit in IPOLZN
c                     if  IPOLZN .NE. 0)

c  NTRM            No. of terms in Mie series

c  PASS1           TRUE on first entry, FALSE thereafter; for self-test

c  PIN(J)          Angular function pi-sub-n ( Ref. 2, Eq. 3 )
c                     at J-th angle
c  PINM1(J)        pi-sub-(n-1) ( see PIn ) at J-th angle

c  PSINM1          Ricatti-Bessel function psi-sub-(n-1), argument XX
c  PSIN            Ricatti-Bessel function psi-sub-n of argument XX
c                     ( Ref. 1, p. 11 ff. )

c  RBIGA(N)        Bessel function ratio A-sub-N (Ref. 2, Eq. 2)
c                     ( REAL version, for when imag refrac index = 0 )

c  RIORIV          1 / Mre

c  RN              1 / N

c  RTMP            (REAL) temporary variable

c  SP(J)           S+  for J-th angle  ( Ref. 2, p. 1507 )
c  SM(J)           S-  for J-TH angle  ( Ref. 2, p. 1507 )
c  SPS(J)          S+  for (NUMANG+1-J)-th angle ( ANYANG=FALSE )
c  SMS(J)          S-  for (NUMANG+1-J)-th angle ( ANYANG=FALSE )

c  TAUN            Angular function tau-sub-n ( Ref. 2, Eq. 4 )
c                     at J-th angle

c  TCOEF           N ( N+1 ) ( 2N+1 ) (for summing TFORW,TBACK series)

c  TWONP1          2N + 1

c  YESANG          TRUE if scattering amplitudes are to be calculated

c  ZETNM1          Ricatti-Bessel function  zeta-sub-(n-1) of argument
c                     XX  ( Ref. 2, Eq. 17 )
c  ZETN            Ricatti-Bessel function  zeta-sub-n of argument XX
c ----------------------------------------------------------------------


      IMPLICIT  NONE

c ----------------------------------------------------------------------
c --------  I / O SPECIFICATIONS FOR SUBROUTINE MIEV0  -----------------
c ----------------------------------------------------------------------
      LOGICAL
     &   ANYANG, PERFCT, PRNT(*)
      INTEGER
     &   IPOLZN, MOMDIM, NUMANG, NMOM
      DOUBLE PRECISION     
     &   GQSC, MIMCUT, PMOM( 0:MOMDIM, * ), QEXT, QSCA, SPIKE,
     &   XMU(*), XX
      DOUBLE COMPLEX  
     &   CREFIN, SFORW, SBACK, S1(*), S2(*), TFORW(*), TBACK(*)
c ----------------------------------------------------------------------

c                                  ** NOTE --  MAXTRM = 10100  is neces-
c                                  ** sary to do some of the test probs,
c                                  ** but 1100 is sufficient for most
c                                  ** conceivable applications
c     .. Parameters ..

      INTEGER
     &   MAXANG, MXANG2
      PARAMETER( 
c     &   MAXANG = 501, 
     &   MAXANG = 180001, 
     &   MXANG2 = MAXANG / 2 + 1 )
      INTEGER
     &   MAXTRM
C      PARAMETER( MAXTRM = 10100 )
      PARAMETER( MAXTRM = 40000 )
      DOUBLE PRECISION
     &   ONE, ONETHR, TWO, ZERO
      PARAMETER ( 
     &   ZERO = 0.D0,
     &   ONETHR = 1.D0 / 3.D0,
     &   ONE = 1.D0,
     &   TWO = 2.D0)

c     ..
c     .. Local Scalars ..

      LOGICAL 
     &   NOABS, PASS1, YESANG
      INTEGER 
     &   I, J, N, NANGD2, NPQUAN, NTRM
      DOUBLE PRECISION      
     &   CHIN, CHINM1, COEFF, DENAN, DENBN, FN, MIM, MM, MRE,
     &   NP1DN, PSIN, PSINM1, RATIO, RIORIV, RN, RTMP, TAUN,
     &   TCOEF, TWONP1, XINV
      DOUBLE COMPLEX   
     &   AN, ANM1, ANP, ANPM, BN, BNM1, BNP, BNPM, CDENAN,
     &   CDENBN, CIOR, CIORIV, CSUM1, CSUM2, CTMP, ZET, 
     &   ZETN, ZETNM1
c     ..
c     .. Local Arrays ..

      LOGICAL 
     &   CALCMO( 4 )
      DOUBLE PRECISION      
     &   PIN( MAXANG ), PINM1( MAXANG ), RBIGA( MAXTRM )
      DOUBLE COMPLEX   
     &   CBIGA( MAXTRM ), LITA( MAXTRM ), LITB( MAXTRM ),
     &   SM( MAXANG ), SMS( MXANG2 ), SP( MAXANG ), SPS( MXANG2 )
c     ..
c     .. External Subroutines ..

      EXTERNAL 
     &   BIGA, CKINMI, ERRMSG, LPCOEF, MIPRNT, SMALL1, SMALL2,
     &   TESTMI

c     .. Intrinsic Functions ..

      INTRINSIC
     &   ABS,CONJG,COS,DCMPLX,DCONJG,DBLE,DIMAG,MAX,MIN,REAL,SIN

      SAVE 
     &   PASS1

c Statement Functions (introduced to enhance portability)

      DOUBLE PRECISION
     &   SQ
      DOUBLE COMPLEX
     &   CJG

      CJG(CTMP)=DCONJG(CTMP)
      SQ(CTMP)=ABS(CTMP)**2

      DATA
     &   PASS1 / .TRUE. /


c                    ** Save some input variables and replace them
c                    ** with values needed to do the self-test

      IF( PASS1 ) CALL TESTMI( .FALSE., XX, CREFIN, MIMCUT, PERFCT,
     &                         ANYANG, NMOM, IPOLZN, NUMANG, XMU, QEXT,
     &                         QSCA, GQSC, SFORW, SBACK, S1, S2, TFORW,
     &                         TBACK, PMOM, MOMDIM )

   10 CONTINUE
c                                        ** Check input and calculate
c                                        ** certain variables from input

      CALL CKINMI( NUMANG, MAXANG, XX, PERFCT, CREFIN, MOMDIM, NMOM,
     &             IPOLZN, ANYANG, XMU, CALCMO, NPQUAN )


      IF( PERFCT .AND. XX.LE.0.1D0 ) THEN
c                                            ** Use totally-reflecting
c                                            ** small-particle limit

         CALL SMALL1( XX, NUMANG, XMU, QEXT, QSCA, GQSC, SFORW, SBACK,
     &                S1, S2, TFORW, TBACK, LITA, LITB )

         NTRM = 2
         GO TO  100

      END IF


      NOABS = .TRUE.

      IF( .NOT.PERFCT ) THEN

         CIOR = CREFIN

         IF( DIMAG(CIOR).GT.ZERO ) CIOR = CJG( CIOR )

         MRE    = DBLE( CIOR )
         MIM    = -DIMAG( CIOR )
         NOABS  = MIM.LE.MIMCUT
         CIORIV = ONE / CIOR
         RIORIV = ONE / MRE

         IF( XX*MAX( ONE, ABS(CIOR) ).LE.0.1D0 ) THEN

c                                    ** Use general-refractive-index
c                                    ** small-particle limit

            CALL SMALL2( XX, CIOR, MIM.GT.MIMCUT, NUMANG, XMU, QEXT,
     &                   QSCA, GQSC, SFORW, SBACK, S1, S2, TFORW,
     &                   TBACK, LITA, LITB )

            NTRM = 2
            GO TO  100

         END IF

      END IF


      NANGD2 = ( NUMANG + 1 ) / 2
      YESANG = NUMANG.GT.0

c                              ** Number of terms in Mie series; Eq R50
      IF( XX.LE.8.D0 ) THEN

         NTRM = INT(XX + 4.0*XX**ONETHR + ONE)

      ELSE IF( XX.LT.4200.D0 ) THEN

         NTRM = INT(XX + 4.05*XX**ONETHR + TWO)

      ELSE

         NTRM = INT(XX + 4.*XX**ONETHR + TWO)

      END IF

      IF( NTRM+1 .GT. MAXTRM )
     &    CALL ERRMSG('MIEV0--PARAMETER MaxTrm TOO SMALL',.TRUE.)

c                            ** Calculate logarithmic derivatives of
c                            ** J-Bessel-fcn., A-sub-(1 to NTrm)

      IF( .NOT.PERFCT ) CALL BIGA( CIOR, XX, NTRM, NOABS, YESANG, RBIGA,
     &                             CBIGA )

c                            ** Initialize Ricatti-Bessel functions
c                            ** (psi,chi,zeta)-sub-(0,1) for upward
c                            ** recurrence ( Eq. R19 )
      XINV   = ONE / XX
      PSINM1 = SIN( XX )
      CHINM1 = COS( XX )
      PSIN   = PSINM1*XINV - CHINM1
      CHIN   = CHINM1*XINV + PSINM1

c 01.11.05 (BTD) replace COMPLEX by DCMPLX
c      ZETNM1 = COMPLEX( PSINM1, CHINM1 )
c      ZETN   = COMPLEX( PSIN, CHIN )

      ZETNM1 = DCMPLX( PSINM1, CHINM1 )
      ZETN   = DCMPLX( PSIN, CHIN )
c                                     ** Initialize previous coeffi-
c                                     ** cients for GQSC series
c 01.11.05 (BTD) replace COMPLEX by DCMPLX
c      ANM1 = COMPLEX( ZERO, ZERO )
c      BNM1 = COMPLEX( ZERO, ZERO )

      ANM1 = DCMPLX( ZERO, ZERO )
      BNM1 = DCMPLX( ZERO, ZERO )
c                             ** Initialize angular function  pi
c                             ** and sums for S+, S- ( Ref. 2, p. 1507 )
      IF( ANYANG ) THEN

         DO 20 J = 1, NUMANG
c                             ** Eq. R39
            PINM1( J ) = ZERO
            PIN( J ) = ONE

c 01.11.05 (BTD) replace COMPLEX by DCMPLX
c            SP( J ) = COMPLEX( ZERO, ZERO )
c            SM( J ) = COMPLEX( ZERO, ZERO )
            SP( J ) = DCMPLX( ZERO, ZERO )
            SM( J ) = DCMPLX( ZERO, ZERO )
   20    CONTINUE

      ELSE

         DO 30 J = 1, NANGD2
c                             ** Eq. R39
            PINM1( J ) = ZERO
            PIN( J ) = ONE

c 01.11.05 (BTD) replace COMPLEX by DCMPLX
c            SP( J ) = COMPLEX( ZERO, ZERO )
c            SM( J ) = COMPLEX( ZERO, ZERO )
c            SPS( J ) = COMPLEX( ZERO, ZERO )
c            SMS( J ) = COMPLEX( ZERO, ZERO )
            SP( J ) = DCMPLX( ZERO, ZERO )
            SM( J ) = DCMPLX( ZERO, ZERO )
            SPS( J ) = DCMPLX( ZERO, ZERO )
            SMS( J ) = DCMPLX( ZERO, ZERO )
   30    CONTINUE

      END IF

c                       ** Initialize Mie sums for efficiencies, etc.
      QSCA  = ZERO
      GQSC  = ZERO

c 01.11.05 (BTD) replace COMPLEX by DCMPLX
c      SFORW = COMPLEX( ZERO, ZERO )
c      SBACK = COMPLEX( ZERO, ZERO )
c      CSUM1 = COMPLEX( ZERO, ZERO )
c      CSUM2 = COMPLEX( ZERO, ZERO )
      SFORW = DCMPLX( ZERO, ZERO )
      SBACK = DCMPLX( ZERO, ZERO )
      CSUM1 = DCMPLX( ZERO, ZERO )
      CSUM2 = DCMPLX( ZERO, ZERO )


c ---------  LOOP TO SUM MIE SERIES  -----------------------------------

      MM     = +ONE
      SPIKE  = ONE

      DO 60  N = 1, NTRM
c                           ** Compute various numerical coefficients
         FN     = N
         RN     = ONE / FN
         NP1DN  = ONE + RN
         TWONP1 = 2*N + 1
         COEFF  = TWONP1 / ( FN * ( N + 1 ) )
         TCOEF  = TWONP1 * ( FN * ( N + 1 ) )

c                           ** Calculate Mie series coefficients
         IF( PERFCT ) THEN
c                                 ** Totally-reflecting case; Eq R/A.1,2

            AN = ( ( FN*XINV )*PSIN - PSINM1 ) /
     &           ( ( FN*XINV )*ZETN - ZETNM1 )
            BN = PSIN / ZETN

         ELSE IF( NOABS ) THEN
c                                      ** No-absorption case; Eq (R16)

            CDENAN = ( RIORIV*RBIGA(N) + ( FN*XINV ) ) * ZETN - ZETNM1
            AN   = ( ( RIORIV*RBIGA(N) + ( FN*XINV ) ) * PSIN - PSINM1 )
     &             / CDENAN
            CDENBN = ( MRE*RBIGA(N) + ( FN*XINV ) ) * ZETN - ZETNM1
            BN   = ( ( MRE*RBIGA(N) + ( FN*XINV ) ) * PSIN - PSINM1 )
     &             / CDENBN

         ELSE
c                                       ** Absorptive case; Eq (R16)

            CDENAN = ( CIORIV*CBIGA( N ) + ( FN*XINV ) )*ZETN - ZETNM1
            CDENBN =   ( CIOR*CBIGA( N ) + ( FN*XINV ) )*ZETN - ZETNM1
            AN   = ( ( CIORIV*CBIGA( N ) + ( FN*XINV ) )*PSIN - PSINM1 )
     &             / CDENAN
            BN     = ( ( CIOR*CBIGA( N ) + ( FN*XINV ) )*PSIN - PSINM1 )
     &             / CDENBN
c                                         ** Eq (R7)

            QSCA   = QSCA + TWONP1*( SQ( AN ) + SQ( BN ) )

         END IF
c                       ** Save Mie coefficients for PMOM calculation

         LITA( N ) = AN
         LITB( N ) = BN


         IF( .NOT.PERFCT .AND. N.GT.XX ) THEN
c                                               ** Flag resonance spikes
            DENAN  = ABS( CDENAN )
            DENBN  = ABS( CDENBN )
c                                                   ** Eq. R/B.9
            RATIO  = DENAN / DENBN
c                                                   ** Eq. R/B.10
            IF( RATIO.LE.0.2 .OR. RATIO.GE.5.0 )
     &          SPIKE = MIN( SPIKE, DENAN, DENBN )

         END IF
c                                  ** Increment Mie sums for non-angle-
c                                  ** dependent quantities

c                                                   ** Eq. R/B.2
         SFORW = SFORW + TWONP1*( AN + BN )
c                                                   ** Eq. R/B.5,6
         CSUM1 = CSUM1 + TCOEF *( AN - BN )
c                                                   ** Eq. R/B.1
         SBACK = SBACK + ( MM*TWONP1 )*( AN - BN )
c                                                   ** Eq. R/B.7,8
         CSUM2 = CSUM2 + ( MM*TCOEF ) *( AN + BN )

c                                         ** Eq (R8)

         GQSC  = GQSC  + ( FN - RN ) * DBLE( ANM1 * CJG( AN ) +
     &                                       BNM1 * CJG( BN ) )
     &           + COEFF * DBLE( AN * CJG( BN ) )


         IF( YESANG ) THEN
c                                      ** Put Mie coefficients in form
c                                      ** needed for computing S+, S-
c                                      ** ( Eq R10 )
            ANP = COEFF*( AN + BN )
            BNP = COEFF*( AN - BN )

c                                      ** Increment Mie sums for S+, S-
c                                      ** while upward recursing
c                                      ** angular functions pi and tau
            IF( ANYANG ) THEN
c                                         ** Arbitrary angles

c                                              ** vectorizable loop
               DO 40 J = 1, NUMANG
c                                                 ** Eq. (R37b)

                  RTMP = ( XMU(J) * PIN(J) ) - PINM1( J )

c                                                 ** Eq. (R38b)
                  TAUN   = FN * RTMP - PINM1( J )

c                                                   ** Eq (R10)

                  SP( J ) = SP( J ) + ANP * ( PIN( J ) + TAUN )
                  SM( J ) = SM( J ) + BNP * ( PIN( J ) - TAUN )

                  PINM1( J ) = PIN( J )
c                                                 ** Eq. R37c

                  PIN( J ) = ( XMU( J ) * PIN( J ) ) + NP1DN * RTMP
   40          CONTINUE

            ELSE
c                                  ** Angles symmetric about 90 degrees
               ANPM = MM*ANP
               BNPM = MM*BNP
c                                          ** vectorizable loop
               DO 50 J = 1, NANGD2
c                                                 ** Eq. (R37b)

                  RTMP = ( XMU(J) * PIN(J) ) - PINM1( J )

c                                                 ** Eq. (R38b)
                  TAUN = FN * RTMP - PINM1( J )

c                                                 ** Eq (R10,12)

                  SP ( J ) = SP ( J ) + ANP * ( PIN( J ) + TAUN )
                  SMS( J ) = SMS( J ) + BNPM *( PIN( J ) + TAUN )
                  SM ( J ) = SM ( J ) + BNP * ( PIN( J ) - TAUN )
                  SPS( J ) = SPS( J ) + ANPM *( PIN( J ) - TAUN )

                  PINM1( J ) = PIN( J )
c                                                 ** Eq. R37c

                  PIN( J ) = ( XMU(J) * PIN(J) ) + NP1DN * RTMP
   50          CONTINUE

            END IF

         END IF
c                          ** Update relevant quantities for next
c                          ** pass through loop
         MM   = - MM
         ANM1 = AN
         BNM1 = BN
c                           ** Upward recurrence for Ricatti-Bessel
c                           ** functions ( Eq. R17 )

         ZET    = ( TWONP1*XINV ) * ZETN - ZETNM1
         ZETNM1 = ZETN
         ZETN   = ZET
         PSINM1 = PSIN
         PSIN   = DBLE( ZETN )

   60 CONTINUE

c ---------- END LOOP TO SUM MIE SERIES --------------------------------


c                                         ** Eq (R6)
      QEXT = TWO / XX**2*DBLE( SFORW )

      IF( PERFCT .OR. NOABS ) THEN

         QSCA = QEXT

      ELSE

         QSCA = TWO/ XX**2 * QSCA

      END IF

      GQSC   = 4./ XX**2 * GQSC
      SFORW  = 0.5*SFORW
      SBACK  = 0.5*SBACK
      TFORW( 1 ) =  0.5*SFORW - 0.125*CSUM1
      TFORW( 2 ) =  0.5*SFORW + 0.125*CSUM1
      TBACK( 1 ) = -0.5*SBACK + 0.125*CSUM2
      TBACK( 2 ) =  0.5*SBACK + 0.125*CSUM2


      IF( YESANG ) THEN
c                                ** Recover scattering amplitudes
c                                ** from S+, S- ( Eq (R11) )

         IF( ANYANG ) THEN
c                                         ** vectorizable loop
            DO 70 J = 1, NUMANG
c                                                   ** Eq (R11)
               S1( J ) = 0.5*( SP( J ) + SM( J ) )
               S2( J ) = 0.5*( SP( J ) - SM( J ) )
   70       CONTINUE

         ELSE
c                                         ** vectorizable loop
            DO 80 J = 1, NANGD2
c                                                   ** Eq (R11)
               S1( J ) = 0.5*( SP( J ) + SM( J ) )
               S2( J ) = 0.5*( SP( J ) - SM( J ) )
   80       CONTINUE
c                                         ** vectorizable loop
            DO 90 J = 1, NANGD2
               S1( NUMANG + 1 - J ) = 0.5*( SPS( J ) + SMS( J ) )
               S2( NUMANG + 1 - J ) = 0.5*( SPS( J ) - SMS( J ) )
   90       CONTINUE

         END IF

      END IF
c                                 ** Calculate Legendre moments

  100 CONTINUE

      IF( NMOM.GT.0 ) CALL LPCOEF( NTRM, NMOM, IPOLZN, MOMDIM, CALCMO,
     &                             NPQUAN, LITA, LITB, PMOM )

      IF( DIMAG( CREFIN ).GT.ZERO ) THEN
c                                         ** Take complex conjugates
c                                         ** of scattering amplitudes

         SFORW = CJG( SFORW )
         SBACK = CJG( SBACK )

         DO 110 I = 1, 2
            TFORW( I ) = CJG( TFORW( I ) )
            TBACK( I ) = CJG( TBACK( I ) )
  110    CONTINUE

         DO 120 J = 1, NUMANG
            S1( J ) = CJG( S1( J ) )
            S2( J ) = CJG( S2( J ) )
  120    CONTINUE

      END IF


      IF( PASS1 ) THEN
c                           ** Compare test case results with
c                           ** correct answers and abort if bad;
c                           ** otherwise restore user input and proceed

         CALL TESTMI( .TRUE., XX, CREFIN, MIMCUT, PERFCT, ANYANG, NMOM,
     &                IPOLZN, NUMANG, XMU, QEXT, QSCA, GQSC, SFORW,
     &                SBACK, S1, S2, TFORW, TBACK, PMOM, MOMDIM )

         PASS1  = .FALSE.
         GO TO  10

      END IF


      IF( PRNT( 1 ) .OR. PRNT( 2 ) ) 
     &  CALL MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT,
     &               QSCA, GQSC, NMOM, IPOLZN, MOMDIM, CALCMO, PMOM,
     &               SFORW, SBACK, TFORW, TBACK, S1, S2 )

      RETURN

      END

      SUBROUTINE CKINMI( NUMANG, MAXANG, XX, PERFCT, CREFIN, MOMDIM,
     &                   NMOM, IPOLZN, ANYANG, XMU, CALCMO, NPQUAN )

c        Check for bad input to MIEV0 and calculate CALCMO, NPQUAN

c     Routines called :  ERRMSG, WRTBAD, WRTDIM


      IMPLICIT NONE

c     .. Parameters ..

      DOUBLE PRECISION ZERO
      PARAMETER(
     &   ZERO = 0.D0 )

c     .. Scalar Arguments ..

      LOGICAL
     &   ANYANG, PERFCT
      INTEGER
     &   IPOLZN, MAXANG, MOMDIM, NMOM, NPQUAN, NUMANG
      DOUBLE PRECISION
     &   XX
      DOUBLE COMPLEX
     &   CREFIN
c     ..
c     .. Array Arguments ..

      LOGICAL 
     &   CALCMO( * )
      DOUBLE PRECISION 
     &   XMU( * )
c     ..
c     .. Local Scalars ..

      CHARACTER 
     &   STRING*4
      LOGICAL
     &   INPERR
      INTEGER
     &   I, IP, J, L
c     ..
c     .. External Functions ..

      LOGICAL
     &   WRTBAD, WRTDIM
      EXTERNAL
     &   WRTBAD, WRTDIM
c     ..
c     .. External Subroutines ..

      EXTERNAL
     &   ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC
     &   ABS, DBLE, DIMAG, ICHAR, MAX
c     ..


      INPERR = .FALSE.

      IF( NUMANG.GT.MAXANG ) INPERR = WRTDIM( 'MaxAng', NUMANG )
      IF( NUMANG.LT.0 ) INPERR = WRTBAD( 'NUMANG' )

      IF( XX.LT.ZERO ) INPERR = WRTBAD( 'XX' )

      IF( .NOT.PERFCT .AND. DBLE( CREFIN ).LE.ZERO )
     &    INPERR = WRTBAD( 'CREFIN' )

      IF( MOMDIM.LT.0 ) INPERR = WRTBAD( 'MOMDIM' )


      IF( NMOM.NE.0 ) THEN

         IF( NMOM.LT.0 .OR. NMOM.GT.MOMDIM ) INPERR = WRTBAD( 'NMOM' )

         IF( ABS( IPOLZN ).GT.4444 ) INPERR = WRTBAD( 'IPOLZN' )

         NPQUAN = 0

         DO 10 L = 1, 4
            CALCMO( L ) = .FALSE.
   10    CONTINUE

         IF( IPOLZN.NE.0 ) THEN
c                                 ** Parse out IPOLZN into its digits
c                                 ** to find which phase quantities are
c                                 ** to have their moments calculated

            WRITE( STRING, '(I4)' ) ABS( IPOLZN )

            DO 20 J = 1, 4
               IP = ICHAR( STRING( J:J ) ) - ICHAR( '0' )

               IF( IP.GE.1 .AND. IP.LE.4 ) CALCMO( IP ) = .TRUE.

               IF( IP.EQ.0 .OR. ( IP.GE.5 .AND. IP.LE.9 ) )
     &             INPERR = WRTBAD( 'IPOLZN' )

               NPQUAN = MAX( NPQUAN, IP )
   20       CONTINUE

         END IF

      END IF


      IF( ANYANG ) THEN
c                                ** Allow for slight imperfections in
c                                ** computation of cosine
         DO 30 I = 1, NUMANG

            IF( XMU( I ).LT.-1.00001 .OR. XMU( I ).GT.1.00001 )
     &          INPERR = WRTBAD( 'XMU' )

   30    CONTINUE

      ELSE

         DO 40 I = 1, ( NUMANG + 1 ) / 2

            IF( XMU( I ).LT.-0.00001 .OR. XMU( I ).GT.1.00001 )
     &          INPERR = WRTBAD( 'XMU' )

   40    CONTINUE

      END IF


      IF( INPERR ) CALL ERRMSG( 'MIEV0--Input error(S).  Aborting...',
     &                          .TRUE. )

c 2002.02.03 (BTD replaced following to get more specific info
c      IF( XX.GT.20000.0 .OR. DBLE( CREFIN ).GT.10.0 .OR.
c     &    ABS( DIMAG( CREFIN ) ).GT.10.0 )
c     &    CALL ERRMSG( 'MIEV0--XX or CREFIN outside tested range',
c     &    .FALSE.)

      IF(XX.GT.20000.0)CALL ERRMSG('MIEV0 -- XX > 20000',.FALSE.)
      IF(DBLE(CREFIN).GT.10.0)CALL ERRMSG('MIEV0 -- Re(m) > 10',.FALSE.)
      IF(DIMAG(CREFIN).GT.10.)
     &   CALL ERRMSG('MIEV0 -- Im(m) > 10',.FALSE.)
      RETURN
      END

      SUBROUTINE LPCOEF( NTRM, NMOM, IPOLZN, MOMDIM, CALCMO, NPQUAN, A,
     &                   B, PMOM )

c         Calculate Legendre polynomial expansion coefficients (also
c         called moments) for phase quantities ( Ref. 5 formulation )

c     INPUT:  NTRM                    Number terms in Mie series
c             NMOM, IPOLZN, MOMDIM    MIEV0 arguments
c             CALCMO                  Flags calculated from IPOLZN
c             NPQUAN                  Defined in MIEV0
c             A, B                    Mie series coefficients

c     OUTPUT: PMOM                   Legendre moments (MIEV0 argument)

c     Routines called :  ERRMSG, LPCO1T, LPCO2T

c     *** NOTES ***

c         (1)  Eqs. 2-5 are in error in Dave, Appl. Opt. 9,
c         1888 (1970).  Eq. 2 refers to M1, not M2;  eq. 3 refers to
c         M2, not M1.  In eqs. 4 and 5, the subscripts on the second
c         term in square brackets should be interchanged.

c         (2)  The general-case logic in this subroutine works correctly
c         in the two-term Mie series case, but subroutine LPCO2T
c         is called instead, for speed.

c         (3)  Subroutine  LPCO1T, to do the one-term case, is never
c         called within the context of MIEV0, but is included for
c         complete generality.

c         (4)  Some improvement in speed is obtainable by combining the
c         310- and 410-loops, if moments for both the third and fourth
c         phase quantities are desired, because the third phase quantity
c         is the real part of a complex series, while the fourth phase
c         quantity is the imaginary part of that very same series.  But
c         most users are not interested in the fourth phase quantity,
c         which is related to circular polarization, so the present
c         scheme is usually more efficient.


c           ** Definitions of local variables ***

c      AM(M)       Numerical coefficients  a-sub-m-super-l
c                     in Dave, Eqs. 1-15, as simplified in Ref. 5.

c      BI(I)       Numerical coefficients  b-sub-i-super-l
c                     in Dave, Eqs. 1-15, as simplified in Ref. 5.

c      BIDEL(I)    1/2 Bi(I) times factor capital-del in Dave

c      CM,DM()     Arrays C and D in Dave, Eqs. 16-17 (Mueller form),
c                     calculated using recurrence derived in Ref. 5

c      CS,DS()     Arrays C and D in Ref. 4, Eqs. A5-A6 (Sekera form),
c                     calculated using recurrence derived in Ref. 5

c      C,D()       Either CM,DM or CS,DS, depending on IPOLZN

c      EVENL       True for even-numbered moments;  false otherwise

c      IDEL        1 + little-del  in Dave

c      MAXTRM      Max. no. of terms in Mie series

c      MAXMOM      Max. no. of non-zero moments

c      NUMMOM      Number of non-zero moments

c      RECIP(K)    1 / K


      IMPLICIT NONE

c     .. Parameters ..

      INTEGER 
     &   MAXTRM, MAXMOM, MXMOM2, MAXRCP
C      PARAMETER(MAXTRM = 1102)
      PARAMETER(MAXTRM = 40002)
      PARAMETER(
     &   MAXMOM = 2*MAXTRM, MXMOM2 = MAXMOM / 2,
     &   MAXRCP = 4*MAXTRM + 2 )
      DOUBLE PRECISION ONE,TWO,ZERO
      PARAMETER(
     &   ZERO = 0.D0,
     &   ONE = 1.D0 ,
     &   TWO = 2.D0 )
      DOUBLE COMPLEX CXZERO
      PARAMETER(
     &   CXZERO = (0.D0,0.D0))
c     ..
c     .. Scalar Arguments ..

      INTEGER 
     &   IPOLZN, MOMDIM, NMOM, NPQUAN, NTRM
c     ..
c     .. Array Arguments ..

      LOGICAL 
     &   CALCMO( * )
      DOUBLE PRECISION 
     &   PMOM( 0:MOMDIM, * )
      DOUBLE COMPLEX 
     &  A( * ), B( * )
c     ..
c     .. Local Scalars ..

      LOGICAL 
     &   EVENL, PASS1
      INTEGER 
     &   I, IDEL, IMAX, J, K, L, LD2, M, MMAX, NUMMOM
      DOUBLE PRECISION 
     &   SUM
      DOUBLE COMPLEX
     &   CTMP
c     ..
c     .. Local Arrays ..

      DOUBLE PRECISION
     &   AM( 0:MAXTRM ), BI( 0:MXMOM2 ), BIDEL( 0:MXMOM2 ),
     &   RECIP( MAXRCP )
      DOUBLE COMPLEX
     &   C( MAXTRM ), CM( MAXTRM ), CS( MAXTRM ), D( MAXTRM ),
     &   DM( MAXTRM ), DS( MAXTRM )
c     ..
c     .. External Subroutines ..

      EQUIVALENCE ( C, CM ), ( D, DM )
c     ..

      SAVE PASS1, RECIP


      EXTERNAL ERRMSG, LPCO1T, LPCO2T

c Intrinsic Functions 

      INTRINSIC CONJG, DCONJG, DIMAG, MAX, MIN, MOD, DBLE

c Statement functions (introduced to enhance portability)

      DOUBLE COMPLEX CJG
      CJG(CTMP)=DCONJG(CTMP)

c     .. Equivalences ..

      DATA PASS1 / .TRUE. /


      IF( PASS1 ) THEN

         DO 10 K = 1, MAXRCP
            RECIP( K ) = ONE / K
   10    CONTINUE

         PASS1  = .FALSE.

      END IF


      DO 30 J = 1, MAX( 1, NPQUAN )

         DO 20 L = 0, NMOM
            PMOM( L, J ) = ZERO
   20    CONTINUE

   30 CONTINUE


      IF( NTRM.EQ.1 ) THEN

         CALL LPCO1T( NMOM, IPOLZN, MOMDIM, CALCMO, A, B, PMOM )

         RETURN

      ELSE IF( NTRM.EQ.2 ) THEN

         CALL LPCO2T( NMOM, IPOLZN, MOMDIM, CALCMO, A, B, PMOM )

         RETURN

      END IF


      IF( NTRM + 2.GT.MAXTRM )
     &    CALL ERRMSG('LPCoef--PARAMETER MaxTrm too small',.TRUE.)

c                                     ** Calculate Mueller C, D arrays
      CM( NTRM + 2 ) = CXZERO
      DM( NTRM + 2 ) = CXZERO
      CM( NTRM + 1 ) = ( ONE - RECIP( NTRM+1 ) ) * B( NTRM )
      DM( NTRM + 1 ) = ( ONE - RECIP( NTRM+1 ) ) * A( NTRM )
      CM( NTRM ) = ( RECIP( NTRM ) + RECIP( NTRM+1 ) ) * A( NTRM ) +
     &             ( ONE - RECIP( NTRM ) )*B( NTRM-1 )
      DM( NTRM ) = ( RECIP( NTRM ) + RECIP( NTRM+1 ) ) * B( NTRM ) +
     &             ( ONE - RECIP( NTRM ) )*A( NTRM-1 )

      DO 40 K = NTRM-1, 2, -1
         CM( K ) = CM( K+2 ) - ( ONE + RECIP(K+1) ) * B( K+1 )
     &                       + ( RECIP(K) + RECIP(K+1) ) * A( K )
     &                       + ( ONE - RECIP(K) ) * B( K-1 )
         DM( K ) = DM( K+2 ) - ( ONE + RECIP(K+1) ) * A( K+1 )
     &                       + ( RECIP(K) + RECIP(K+1) ) * B( K )
     &                       + ( ONE - RECIP(K) ) * A( K-1 )
   40 CONTINUE

      CM( 1 ) = CM( 3 ) + 1.5D0 * ( A( 1 ) - B( 2 ) )
      DM( 1 ) = DM( 3 ) + 1.5D0 * ( B( 1 ) - A( 2 ) )


      IF( IPOLZN.GE.0 ) THEN

         DO 50 K = 1, NTRM + 2
            C( K ) = ( 2*K - 1 ) * CM( K )
            D( K ) = ( 2*K - 1 ) * DM( K )
   50    CONTINUE

      ELSE
c                                    ** Compute Sekera C and D arrays
         CS( NTRM + 2 ) = CXZERO
         DS( NTRM + 2 ) = CXZERO
         CS( NTRM + 1 ) = CXZERO
         DS( NTRM + 1 ) = CXZERO

         DO 60 K = NTRM, 1, -1
            CS( K ) = CS( K+2 ) + ( 2*K + 1 ) * ( CM( K+1 ) - B( K ) )
            DS( K ) = DS( K+2 ) + ( 2*K + 1 ) * ( DM( K+1 ) - A( K ) )
   60    CONTINUE

         DO 70 K = 1, NTRM + 2
            C( K ) = ( 2*K - 1 ) * CS( K )
            D( K ) = ( 2*K - 1 ) * DS( K )
   70    CONTINUE

      END IF


      IF( IPOLZN.LT.0 ) NUMMOM = MIN( NMOM, 2*NTRM - 2 )
      IF( IPOLZN.GE.0 ) NUMMOM = MIN( NMOM, 2*NTRM )

      IF( NUMMOM.GT.MAXMOM )
     &    CALL ERRMSG('LPCoef--PARAMETER MaxTrm too small',.TRUE.)


c                          ** Loop over moments

      DO 240 L = 0, NUMMOM

         LD2 = L / 2
         EVENL  = MOD( L, 2 ).EQ.0
c                                    ** Calculate numerical coefficients
c                                    ** a-sub-m and b-sub-i in Dave
c                                    ** double-sums for moments
         IF( L.EQ.0 ) THEN

            IDEL = 1

            DO 80 M = 0, NTRM
               AM( M ) = TWO * RECIP( 2*M + 1 )
   80       CONTINUE

            BI( 0 ) = ONE

         ELSE IF( EVENL ) THEN

            IDEL = 1

            DO 90 M = LD2, NTRM
               AM( M ) = ( ONE + RECIP( 2*M - L + 1 ) ) * AM( M )
   90       CONTINUE

            DO 100 I = 0, LD2 - 1
               BI( I ) = ( ONE - RECIP( L - 2*I ) ) * BI( I )
  100       CONTINUE

            BI( LD2 ) = ( TWO - RECIP( L ) ) * BI( LD2 - 1 )

         ELSE

            IDEL = 2

            DO 110 M = LD2, NTRM
               AM( M ) = ( ONE - RECIP( 2*M + L + 2 ) ) * AM( M )
  110       CONTINUE

            DO 120 I = 0, LD2
               BI( I ) = ( ONE - RECIP( L + 2*I + 1 ) ) * BI( I )
  120       CONTINUE

         END IF
c                                     ** Establish upper limits for sums
c                                     ** and incorporate factor capital-
c                                     ** del into b-sub-i
         MMAX = NTRM - IDEL
         IF( IPOLZN.GE.0 ) MMAX = MMAX + 1
         IMAX = MIN( LD2, MMAX - LD2 )

         IF( IMAX.LT.0 ) GO TO  250

         DO 130 I = 0, IMAX
            BIDEL( I ) = BI( I )
  130    CONTINUE

         IF( EVENL ) BIDEL( 0 ) = 0.5*BIDEL( 0 )

c                                    ** Perform double sums just for
c                                    ** phase quantities desired by user
         IF( IPOLZN.EQ.0 ) THEN

            DO 150 I = 0, IMAX
c                                           ** vectorizable loop

               SUM = ZERO

               DO 140 M = LD2, MMAX - I
                  SUM = SUM + AM( M ) *
     &                      ( DBLE( C(M-I+1) * CJG( C(M+I+IDEL) ) )
     &                      + DBLE( D(M-I+1) * CJG( D(M+I+IDEL) ) ) )
  140          CONTINUE

               PMOM( L, 1 ) = PMOM( L, 1 ) + BIDEL( I ) * SUM

  150       CONTINUE

            PMOM( L, 1 ) = 0.5D0*PMOM( L, 1 )
            GO TO  240

         END IF


         IF( CALCMO( 1 ) ) THEN

            DO 170 I = 0, IMAX

               SUM = ZERO
c                                           ** vectorizable loop
               DO 160 M = LD2, MMAX - I
                  SUM = SUM + AM( M ) *
     &                        DBLE( C(M-I+1) * CJG( C(M+I+IDEL) ) )
  160          CONTINUE

               PMOM( L, 1 ) = PMOM( L, 1 ) + BIDEL( I ) * SUM

  170       CONTINUE

         END IF


         IF( CALCMO( 2 ) ) THEN

            DO 190 I = 0, IMAX

               SUM = ZERO
c                                           ** vectorizable loop
               DO 180 M = LD2, MMAX - I
                  SUM = SUM + AM( M ) *
     &                        DBLE( D(M-I+1) * CJG( D(M+I+IDEL) ) )
  180          CONTINUE

               PMOM( L, 2 ) = PMOM( L, 2 ) + BIDEL( I ) * SUM

  190       CONTINUE

         END IF


         IF( CALCMO( 3 ) ) THEN

            DO 210 I = 0, IMAX

               SUM = ZERO
c                                           ** vectorizable loop
               DO 200 M = LD2, MMAX - I
                  SUM = SUM + AM( M ) *
     &                      ( DBLE( C(M-I+1) * CJG( D(M+I+IDEL) ) )
     &                      + DBLE( C(M+I+IDEL) * CJG( D(M-I+1) ) ) )
  200          CONTINUE

               PMOM( L, 3 ) = PMOM( L, 3 ) + BIDEL( I ) * SUM

  210       CONTINUE

            PMOM( L, 3 ) = 0.5D0*PMOM( L, 3 )
         END IF


         IF( CALCMO( 4 ) ) THEN

            DO 230 I = 0, IMAX

               SUM= ZERO
c                                           ** vectorizable loop
               DO 220 M = LD2, MMAX - I
                  SUM = SUM + AM( M ) *
     &                      ( DIMAG( C(M-I+1) * CJG( D(M+I+IDEL) ) )
     &                      + DIMAG( C(M+I+IDEL) * CJG( D(M-I+1) ) ) )
  220          CONTINUE

               PMOM( L, 4 ) = PMOM( L, 4 ) + BIDEL( I ) * SUM

  230       CONTINUE

            PMOM( L, 4 ) = - 0.5D0 * PMOM( L, 4 )

         END IF

  240 CONTINUE


  250 CONTINUE

      RETURN
      END

      SUBROUTINE LPCO1T( NMOM, IPOLZN, MOMDIM, CALCMO, A, B, PMOM )

c         Calculate Legendre polynomial expansion coefficients (also
c         called moments) for phase quantities in special case where
c         no. terms in Mie series = 1

c        INPUT:  NMOM, IPOLZN, MOMDIM     MIEV0 arguments
c                CALCMO                   Flags calculated from IPOLZN
c                A(1), B(1)               Mie series coefficients

c        OUTPUT: PMOM                     Legendre moments


      IMPLICIT  NONE

c     .. Parameters

      DOUBLE PRECISION ZERO
      PARAMETER(
     &   ZERO = 0.D0 )

c     .. Scalar Arguments ..

      INTEGER
     &   IPOLZN, MOMDIM, NMOM
c     ..
c     .. Array Arguments ..

      LOGICAL
     &   CALCMO( * )
      DOUBLE PRECISION
     &   PMOM( 0:MOMDIM, * )
      DOUBLE COMPLEX
     &   A( * ), B( * )
c     ..
c     .. Local Scalars ..

      INTEGER
     &   L, NUMMOM
      DOUBLE PRECISION
     &   A1SQ, B1SQ
      DOUBLE COMPLEX
     &   A1B1C, CTMP
c     ..
c     .. Intrinsic Functions ..

c 01.11.05 (BTD) attempt to make portable
c      INTRINSIC DIMAG, CONJG, MIN, DBLE, REALPART
      INTRINSIC DIMAG, CONJG, DCONJG, MIN, DBLE
c     ..
c     .. Statement Functions ..

      DOUBLE PRECISION SQ
      DOUBLE COMPLEX CJG
c     ..
c     .. Statement Function definitions ..

c 01.11.05 (BTD) attempt to make portable
c g77 is unhappy with following statement
c      SQ( CTMP ) = DBLE( CTMP )**2 + DIMAG( CTMP )**2
c solaris f77 is unhappy with following statement
c      SQ( CTMP ) = REALPART( CTMP )**2 + DIMAG( CTMP )**2
c try following:

      CJG(CTMP)=DCONJG(CTMP)
      SQ(CTMP)=ABS(CTMP)**2

      A1SQ   = SQ( A( 1 ) )
      B1SQ   = SQ( B( 1 ) )
      A1B1C  = A( 1 ) * CJG( B( 1 ) )


      IF( IPOLZN.LT.0 ) THEN

         IF( CALCMO( 1 ) ) PMOM( 0, 1 ) = 2.25D0*B1SQ

         IF( CALCMO( 2 ) ) PMOM( 0, 2 ) = 2.25D0*A1SQ

         IF( CALCMO( 3 ) ) PMOM( 0, 3 ) = 2.25D0*DBLE( A1B1C )

         IF( CALCMO( 4 ) ) PMOM( 0, 4 ) = 2.25D0*DIMAG( A1B1C )

      ELSE

         NUMMOM = MIN( NMOM, 2 )

c                             ** Loop over moments
         DO 10  L = 0, NUMMOM

            IF( IPOLZN.EQ.0 ) THEN

               IF( L.EQ.0 ) PMOM( L, 1 ) = 1.5D0*( A1SQ + B1SQ )

               IF( L.EQ.1 ) PMOM( L, 1 ) = 1.5D0*DBLE( A1B1C )

               IF( L.EQ.2 ) PMOM( L, 1 ) = 0.15D0*( A1SQ + B1SQ )

               GO TO  10

            END IF


            IF( CALCMO( 1 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 1 ) = 2.25D0*( A1SQ + B1SQ / 3.D0)

               IF( L.EQ.1 ) PMOM( L, 1 ) = 1.5D0*DBLE( A1B1C )

               IF( L.EQ.2 ) PMOM( L, 1 ) = 0.3D0*B1SQ

            END IF


            IF( CALCMO( 2 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 2 ) = 2.25D0*( B1SQ + A1SQ / 3.D0 )

               IF( L.EQ.1 ) PMOM( L, 2 ) = 1.5D0*DBLE( A1B1C )

               IF( L.EQ.2 ) PMOM( L, 2 ) = 0.3D0*A1SQ

            END IF


            IF( CALCMO( 3 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 3 ) = 3.0D0*DBLE( A1B1C )

               IF( L.EQ.1 ) PMOM( L, 3 ) = 0.75D0*( A1SQ + B1SQ )

               IF( L.EQ.2 ) PMOM( L, 3 ) = 0.3D0*DBLE( A1B1C )

            END IF


            IF( CALCMO( 4 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 4 ) = -1.5D0*DIMAG( A1B1C )

               IF( L.EQ.1 ) PMOM( L, 4 ) = ZERO

               IF( L.EQ.2 ) PMOM( L, 4 ) = 0.3D0*DIMAG( A1B1C )

            END IF


   10    CONTINUE

      END IF

      RETURN
      END

      SUBROUTINE LPCO2T( NMOM, IPOLZN, MOMDIM, CALCMO, A, B, PMOM )

c         Calculate Legendre polynomial expansion coefficients (also
c         called moments) for phase quantities in special case where
c         no. terms in Mie series = 2

c        INPUT:  NMOM, IPOLZN, MOMDIM     MIEV0 arguments
c                CALCMO                   Flags calculated from IPOLZN
c                A(1-2), B(1-2)           Mie series coefficients

c        OUTPUT: PMOM                     Legendre moments


      IMPLICIT  NONE

c     .. Parameters ..

      DOUBLE PRECISION FIVTHRD,TENTHRD,TWOTHRD
      PARAMETER(
     &   TWOTHRD = 2.D0/3.D0,
     &   FIVTHRD = 5.D0/3.D0,
     &   TENTHRD = 10.D0/3.D0)
     

c     .. Scalar Arguments ..

      INTEGER   IPOLZN, MOMDIM, NMOM
c     ..
c     .. Array Arguments ..

      LOGICAL   CALCMO( * )
      DOUBLE PRECISION      PMOM( 0:MOMDIM, * )
      DOUBLE COMPLEX   A( * ), B( * )
c     ..
c     .. Local Scalars ..

      INTEGER 
     &   L, NUMMOM
      DOUBLE PRECISION
     &   A2SQ, B2SQ, PM1, PM2
      DOUBLE COMPLEX 
     &   A2C, B2C, CA, CAC, CAT, CB, CBC, CBT, CG, CH, CTMP
c     ..
c     .. Intrinsic Functions ..


c 01.11.05 (BTD) attempt to make portable:
c     &   DIMAG, CONJG, MIN, DBLE, REALPART

      INTRINSIC
     &   DIMAG, DCONJG, MIN, DBLE
c     ..
c     .. Statement Functions ..

      DOUBLE PRECISION SQ
      DOUBLE COMPLEX CJG
c     ..
c     .. Statement Function definitions ..

C 01.11.05 (BTD) attempt to make portable
c g77 is unhappy with following statement:
c      SQ( CTMP ) = DBLE( CTMP )**2 + DIMAG( CTMP )**2
c solaris f77 is unhappy with following:
c      SQ( CTMP ) = REALPART( CTMP )**2 + DIMAG( CTMP )**2
c try following:

      CJG(CTMP)=DCONJG(CTMP)
      SQ(CTMP)=ABS(CTMP)**2

      CA   = 3.D0*A( 1 ) - 5.D0*B( 2 )
      CAT  = 3.D0*B( 1 ) - 5.D0*A( 2 )
      CAC  = CJG( CA )
      A2SQ = SQ( A( 2 ) )
      B2SQ = SQ( B( 2 ) )
      A2C  = CJG( A( 2 ) )
      B2C  = CJG( B( 2 ) )


      IF( IPOLZN.LT.0 ) THEN

c                                   ** Loop over Sekera moments
         NUMMOM = MIN( NMOM, 2 )

         DO 10 L = 0, NUMMOM

            IF( CALCMO( 1 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 1 ) = 0.25D0 * ( SQ( CAT )
     &                                      + (100.D0/3.D0)* B2SQ )

               IF( L.EQ.1 ) PMOM( L, 1 ) = FIVTHRD*DBLE( CAT*B2C )

               IF( L.EQ.2 ) PMOM( L, 1 ) = TENTHRD*B2SQ

            END IF


            IF( CALCMO( 2 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 2 ) = 0.25D0 * ( SQ( CA )
     &                                      + (100.D0/3.D0) * A2SQ )

               IF( L.EQ.1 ) PMOM( L, 2 ) = FIVTHRD*DBLE( CA*A2C )

               IF( L.EQ.2 ) PMOM( L, 2 ) = TENTHRD*A2SQ

            END IF


            IF( CALCMO( 3 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 3 ) = 0.25D0 * DBLE( CAT * CAC
     &                                    + (100.D0/3.D0) * B(2) * A2C )

               IF( L.EQ.1 ) PMOM( L, 3 ) = 5.D0/6.D0*
     &                                     DBLE( B(2)*CAC + CAT*A2C )

               IF( L.EQ.2 ) PMOM( L, 3 ) = TENTHRD* DBLE( B(2)*A2C )

            END IF


            IF( CALCMO( 4 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 4 ) = -0.25D0 * DIMAG( CAT * CAC
     &                                    + (100.D0/3.D0)* B(2) * A2C )

               IF( L.EQ.1 ) PMOM( L, 4 ) = -(5.D0/ 6.D0)*
     &                                     DIMAG( B(2)*CAC + CAT*A2C )

               IF( L.EQ.2 ) PMOM( L, 4 ) = -TENTHRD * DIMAG( B(2)*A2C )

            END IF


   10    CONTINUE


      ELSE

         CB  = 3.D0*B( 1 ) + 5.D0*A( 2 )
         CBT = 3.D0*A( 1 ) + 5.D0*B( 2 )
         CBC = CJG( CB )
         CG  = ( CBC*CBT + 10.D0*( CAC*A( 2 ) + B2C*CAT ) ) / 3.D0
         CH  = 2.D0*( CBC*A( 2 ) + B2C*CBT )

c                               ** Loop over Mueller moments
         NUMMOM = MIN( NMOM, 4 )

         DO 20 L = 0, NUMMOM


            IF( IPOLZN.EQ.0 .OR. CALCMO( 1 ) ) THEN

               IF( L.EQ.0 ) PM1 = 0.25D0*SQ( CA ) + SQ( CB ) / 12.D0
     &                            + FIVTHRD*DBLE( CA*B2C ) + 5.D0*B2SQ

               IF( L.EQ.1 ) PM1 = DBLE( CB * ( CAC / 6.D0+ B2C ) )

               IF( L.EQ.2 ) PM1 = SQ( CB ) / 30.D0+ (20.D0/7.D0)*B2SQ
     &                            + TWOTHRD*DBLE( CA*B2C )

               IF( L.EQ.3 ) PM1 = (2.D0/7.D0) * DBLE( CB*B2C )

               IF( L.EQ.4 ) PM1 = (40.D0/63.D0) * B2SQ

               IF( CALCMO( 1 ) ) PMOM( L, 1 ) = PM1

            END IF


            IF( IPOLZN.EQ.0 .OR. CALCMO( 2 ) ) THEN

               IF( L.EQ.0 ) PM2 = 0.25D0*SQ( CAT ) + SQ( CBT ) / 12.D0
     &                           + FIVTHRD * DBLE( CAT*A2C )
     &                           + 5.D0*A2SQ

               IF( L.EQ.1 ) PM2 = DBLE( CBT *
     &                                 ( CJG( CAT ) / 6.D0+ A2C ) )

               IF( L.EQ.2 ) PM2 = SQ( CBT ) / 30.D0
     &                            + ( 20.D0/7.D0) * A2SQ
     &                            + TWOTHRD * DBLE( CAT*A2C )

               IF( L.EQ.3 ) PM2 = (2.D0/7.D0) * DBLE( CBT*A2C )

               IF( L.EQ.4 ) PM2 = (40.D0/63.D0) * A2SQ

               IF( CALCMO( 2 ) ) PMOM( L, 2 ) = PM2

            END IF


            IF( IPOLZN.EQ.0 ) THEN

               PMOM( L, 1 ) = 0.5D0*( PM1 + PM2 )
               GO TO  20

            END IF


            IF( CALCMO( 3 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 3 ) = 0.25D0 * DBLE( CAC*CAT + CG
     &                                         + 20.D0* B2C * A(2) )

               IF( L.EQ.1 ) PMOM( L, 3 ) = DBLE( CAC*CBT + CBC*CAT
     &                                          + 3.D0*CH ) / 12.D0

               IF( L.EQ.2 ) PMOM( L, 3 ) = 0.1D0 * DBLE( CG
     &                                    + (200.D0/7.D0) * B2C * A(2) )

               IF( L.EQ.3 ) PMOM( L, 3 ) = DBLE( CH ) / 14.D0

               IF( L.EQ.4 ) PMOM( L, 3 ) = (40.D0/63.D0)*DBLE(B2C*A(2) )

            END IF


            IF( CALCMO( 4 ) ) THEN

               IF( L.EQ.0 ) PMOM( L, 4 ) = 0.25D0 * DIMAG( CAC*CAT + CG
     &                                      + 20.D0* B2C * A(2) )

               IF( L.EQ.1 ) PMOM( L, 4 ) = DIMAG( CAC*CBT + CBC*CAT
     &                                          + 3.D0*CH ) / 12.D0

               IF( L.EQ.2 ) PMOM( L, 4 ) = 0.1D0 * DIMAG( CG
     &                                   + (200.D0/7.D0) * B2C * A(2) )

               IF( L.EQ.3 ) PMOM( L, 4 ) = DIMAG( CH ) / 14.D0

               IF( L.EQ.4 ) PMOM( L, 4 ) = 40.D0/63.D0* DIMAG(B2C*A(2))

            END IF


   20    CONTINUE

      END IF

      RETURN
      END

      SUBROUTINE BIGA( CIOR, XX, NTRM, NOABS, YESANG, RBIGA, CBIGA )

c        Calculate logarithmic derivatives of J-Bessel-function

c     Input :  CIOR, XX, NTRM, NOABS, YESANG  (defined in MIEV0)

c    Output :  RBIGA or CBIGA  (defined in MIEV0)

c    Routines called :  CONFRA


c    INTERNAL VARIABLES :

c       CONFRA     Value of Lentz continued fraction for CBIGA(NTRM),
c                     used to initialize downward recurrence

c       DOWN       = True, use down-recurrence.  False, do not.

c       F1,F2,F3   Arithmetic statement functions used in determining
c                     whether to use up-  or down-recurrence
c                     ( Ref. 2, Eqs. 6-8 )

c       MRE        Real refractive index
c       MIM        Imaginary refractive index

c       REZINV     1 / ( MRE * XX ); temporary variable for recurrence
c       ZINV       1 / ( CIOR * XX ); temporary variable for recurrence


      IMPLICIT  NONE

c     .. Parameters ..

      DOUBLE PRECISION ONE
      PARAMETER(
     &   ONE=1.D0)

c      DOUBLE PRECISION ONE,TWO,ZERO
c     &   ZERO = 0.D0,
c     &   ONE = 1.D0,
c     &   TWO = 2.D0)
      DOUBLE COMPLEX CXONEI,CXTWOI
      PARAMETER(
     &   CXONEI=(0.D0,1.D0),
     &   CXTWOI=(0.D0,2.D0))

c     .. Scalar Arguments ..

      LOGICAL   NOABS, YESANG
      INTEGER   NTRM
      DOUBLE PRECISION      XX
      DOUBLE COMPLEX   CIOR
c     ..
c     .. Array Arguments ..

      DOUBLE PRECISION      RBIGA( * )
      DOUBLE COMPLEX   CBIGA( * )
c     ..
c     .. Local Scalars ..

      LOGICAL   DOWN
      INTEGER   N
      DOUBLE PRECISION      MIM, MRE, REZINV, RTMP
      DOUBLE COMPLEX   CTMP, ZINV
c     ..
c     .. External Functions ..

      DOUBLE COMPLEX   CONFRA
      EXTERNAL  CONFRA
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, DIMAG, COS, EXP, DBLE, SIN
c     ..
c     .. Statement Functions ..

c      DOUBLE PRECISION      F1, F2, F3
      DOUBLE PRECISION F1,F2
c     ..
c     .. Statement Function definitions ..

c                                                   ** Eq. R47c
      F1( MRE ) = -8.0 + MRE**2*( 26.22 +
     &            MRE*( -0.4474 + MRE**3*( 0.00204 - 0.000175*MRE ) ) )

c                                                   ** Eq. R47b
      F2( MRE ) = 3.9 + MRE*( -10.8 + 13.78*MRE )
c                                                   ** Eq. R47a
c      F3( MRE ) = -15.04 + MRE*( 8.42 + 16.35*MRE )
c     ..

c                                  ** Decide whether BigA can be
c                                  ** calculated by up-recurrence
      MRE = DBLE( CIOR )
      MIM = ABS( DIMAG( CIOR ) )

      IF( MRE.LT.ONE .OR. MRE.GT.10.0 .OR. MIM.GT.10.0 ) THEN

         DOWN = .TRUE.

      ELSE IF( YESANG ) THEN

         DOWN = .TRUE.
c                                                    ** Eq. R48
         IF( MIM*XX .LT. F2( MRE ) ) DOWN = .FALSE.

      ELSE

         DOWN = .TRUE.
c                                                    ** Eq. R48
         IF( MIM*XX .LT. F1( MRE ) ) DOWN = .FALSE.

      END IF


      ZINV   = ONE / ( CIOR*XX )
      REZINV = ONE / ( MRE*XX )


      IF( DOWN ) THEN
c                          ** Compute initial high-order BigA using
c                          ** Lentz method ( Ref. 1, pp. 17-20 )

         CTMP = CONFRA( NTRM, ZINV )

c                                   *** Downward recurrence for BigA
         IF( NOABS ) THEN
c                                        ** No-absorption case; Eq (R23)
            RBIGA( NTRM ) = DBLE( CTMP )

            DO 10 N = NTRM, 2, -1
               RBIGA( N - 1 ) = ( N*REZINV ) -
     &                          ONE / ( ( N*REZINV ) + RBIGA( N ) )
   10       CONTINUE

         ELSE
c                                         ** Absorptive case; Eq (R23)
            CBIGA( NTRM ) = CTMP

            DO 20 N = NTRM, 2, -1
               CBIGA( N-1 ) = (N*ZINV) - ONE / ( (N*ZINV) + CBIGA( N ) )
   20       CONTINUE

         END IF


      ELSE
c                            *** Upward recurrence for BigA
         IF( NOABS ) THEN
c                                  ** No-absorption case; Eq (R20,21)
            RTMP = SIN( MRE*XX )
            RBIGA( 1 ) = - REZINV + RTMP /
     &                   ( RTMP*REZINV - COS( MRE*XX ) )

            DO 30 N = 2, NTRM
               RBIGA( N ) = -( N*REZINV ) +
     &                      ONE / ( ( N*REZINV ) - RBIGA( N - 1 ) )
   30       CONTINUE

         ELSE
c                                     ** Absorptive case; Eq (R20,22)

            CTMP = EXP( - CXTWOI*CIOR*XX )
            CBIGA( 1 ) = - ZINV + (ONE-CTMP) /
     &             ( ZINV * (ONE-CTMP) - CXONEI*(ONE+CTMP) )

            DO 40 N = 2, NTRM
               CBIGA( N ) = - (N*ZINV) + ONE / ((N*ZINV) - CBIGA( N-1 ))
   40       CONTINUE

         END IF

      END IF

      RETURN
      END

      DOUBLE COMPLEX FUNCTION CONFRA( N, ZINV )

c         Compute Bessel function ratio A-sub-N from its
c         continued fraction using Lentz method

c         ZINV = Reciprocal of argument of A


c    I N T E R N A L    V A R I A B L E S
c    ------------------------------------

c    CAK      Term in continued fraction expansion of A (Eq. R25)

c    CAPT     Factor used in Lentz iteration for A (Eq. R27)

c    CNUMER   Numerator   in capT  ( Eq. R28A )
c    CDENOM   Denominator in capT  ( Eq. R28B )

c    CDTD     Product of two successive denominators of capT factors
c                 ( Eq. R34C )
c    CNTN     Product of two successive numerators of capT factors
c                 ( Eq. R34B )

c    EPS1     Ill-conditioning criterion
c    EPS2     Convergence criterion

c    KK       Subscript k of cAk  ( Eq. R25B )

c    KOUNT    Iteration counter ( used to prevent infinite looping )

c    MAXIT    Max. allowed no. of iterations

c    MM       + 1  and - 1, alternately
c --------------------------------------------------------------------

      IMPLICIT  NONE

c     .. Parameters ..

      DOUBLE PRECISION ONE
      PARAMETER(
     &   ONE=1.D0)

c     .. Scalar Arguments ..

      INTEGER   N
      DOUBLE COMPLEX   ZINV
c     ..
c     .. Local Scalars ..

      INTEGER   KK, KOUNT, MAXIT, MM
      DOUBLE PRECISION      EPS1, EPS2
      DOUBLE COMPLEX   CAK, CAPT, CDENOM, CDTD, CNTN, CNUMER
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, DIMAG, DBLE
c     ..
      DATA      EPS1 / 1.E-2 / , EPS2 / 1.E-8 /
      DATA      MAXIT / 10000 /


c                                 ** Eq. R25a
      CONFRA = ( N + 1 ) * ZINV
      MM     = - 1
      KK     = 2*N + 3
c                                 ** Eq. R25b, k=2
      CAK    = ( MM*KK ) * ZINV
      CDENOM = CAK
      CNUMER = CDENOM + ONE / CONFRA
      KOUNT  = 1

   10 CONTINUE
      KOUNT = KOUNT + 1

      IF( KOUNT.GT.MAXIT )
     &    CALL ERRMSG('ConFra--Iteration failed to converge',.TRUE.)

      MM  = - MM
      KK  = KK + 2
c                                 ** Eq. R25b
      CAK = ( MM*KK ) * ZINV
c                                          ** Eq. R32
      IF( ABS( CNUMER / CAK ).LE.EPS1 .OR.
     &    ABS( CDENOM / CAK ).LE.EPS1 ) THEN

c                                  ** Ill-conditioned case -- stride
c                                  ** two terms instead of one

c                                       ** Eq. R34
         CNTN   = CAK * CNUMER + ONE
         CDTD   = CAK * CDENOM + ONE
c                                           ** Eq. R33
         CONFRA = ( CNTN / CDTD ) * CONFRA

         MM  = - MM
         KK  = KK + 2
c                                 ** Eq. R25b
         CAK = ( MM*KK ) * ZINV
c                                      ** Eq. R35
         CNUMER = CAK + CNUMER / CNTN
         CDENOM = CAK + CDENOM / CDTD
         KOUNT  = KOUNT + 1
         GO TO  10

      ELSE
c                           *** Well-conditioned case

c                                  ** Eq. R27
         CAPT   = CNUMER / CDENOM
c                                  ** Eq. R26
         CONFRA = CAPT * CONFRA
c                                  ** Check for convergence; Eq. R31

         IF (      ABS( DBLE (CAPT) - ONE ).GE.EPS2
     &        .OR. ABS( DIMAG(CAPT) )      .GE.EPS2 )  THEN

c                                        ** Eq. R30
            CNUMER = CAK + ONE / CNUMER
            CDENOM = CAK + ONE / CDENOM

            GO TO  10

         END IF

      END IF


      RETURN

      END

      SUBROUTINE MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT,
     &                   QSCA, GQSC, NMOM, IPOLZN, MOMDIM, CALCMO, PMOM,
     &                   SFORW, SBACK, TFORW, TBACK, S1, S2 )

c         Print scattering quantities of a single particle


      IMPLICIT  NONE

c     .. Scalar Arguments ..

      LOGICAL
     &   PERFCT
      INTEGER
     &   IPOLZN, MOMDIM, NMOM, NUMANG
      DOUBLE PRECISION
     &   GQSC, QEXT, QSCA, XX
      DOUBLE COMPLEX
     &   CREFIN, CTMP, SBACK, SFORW
c     ..
c     .. Array Arguments ..

      LOGICAL
     &   CALCMO( * ), PRNT( * )
      DOUBLE PRECISION
     &   PMOM( 0:MOMDIM, * ), XMU( * )
      DOUBLE COMPLEX
     &   S1( * ), S2( * ), TBACK( * ), TFORW( * )
c     ..
c     .. Local Scalars ..

      CHARACTER FMAT*22
      INTEGER   I, J, M
      DOUBLE PRECISION      FNORM, I1, I2
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC DIMAG, DCONJG, DBLE
c     ..
c statement functions

      DOUBLE COMPLEX CJG
      CJG(CTMP)=DCONJG(CTMP)


      IF( PERFCT ) WRITE( *, '(''1'',10X,A,1P,E11.4)' )
     &    'Perfectly Conducting Case, size parameter =', XX

      IF( .NOT.PERFCT ) WRITE( *, '(''1'',10X,3(A,1P,E11.4))' )
     &    'Refractive Index:  Real ', DBLE( CREFIN ), '  Imag ',
     &    DIMAG( CREFIN ), ',   Size Parameter =', XX


      IF( PRNT( 1 ) .AND. NUMANG.GT.0 ) THEN

         WRITE( *, '(/,A)' )
     &      '    cos(angle)  ------- S1 ---------  ------- S2 ---------'
     &      // '  --- S1*conjg(S2) ---   i1=S1**2   i2=S2**2  (i1+i2)/2'
     &      // '  DEG POLZN'

         DO 10 I = 1, NUMANG
            I1 = DBLE( S1( I ) )**2 + DIMAG( S1( I ) )**2
            I2 = DBLE( S2( I ) )**2 + DIMAG( S2( I ) )**2
            WRITE( *, '( I4, F10.6, 1P,10E11.3 )'   )
     &              I, XMU(I), S1(I), S2(I), S1(I)*CJG(S2(I)),
     &              I1, I2, 0.5*(I1+I2), (I2-I1)/(I2+I1)
   10    CONTINUE

      END IF


      IF( PRNT( 2 ) ) THEN

         WRITE ( *, '(/,A,9X,A,17X,A,17X,A,/,(0P,F7.2, 1P,6E12.3) )' )
     &           '  Angle', 'S-sub-1', 'T-sub-1', 'T-sub-2',
     &               0.0,     SFORW,    TFORW(1),  TFORW(2),
     &              180.,     SBACK,    TBACK(1),  TBACK(2)
         WRITE ( *, '(/,4(A,1P,E11.4))' )
     &           ' Efficiency Factors,  extinction:', QEXT,
     &                              '   scattering:', QSCA,
     &                              '   absorption:', QEXT-QSCA,
     &                           '   rad. pressure:', QEXT-GQSC

         IF( NMOM.GT.0 ) THEN

            WRITE( *, '(/,A)' ) ' Normalized moments of :'

            IF( IPOLZN.EQ.0 ) WRITE( *, '(''+'',27X,A)' )
     &          'Phase Fcn'

            IF( IPOLZN.GT.0 ) WRITE( *, '(''+'',33X,A)' )
     &          'M1           M2          S21          D21'

            IF( IPOLZN.LT.0 ) WRITE( *, '(''+'',33X,A)' )
     &          'R1           R2           R3           R4'

            FNORM = 4./ ( XX**2 * QSCA )

            DO 30  M = 0, NMOM

               WRITE( *, '(A,I4)' ) '      Moment no.', M

               DO 20 J = 1, 4

                  IF( CALCMO( J ) ) THEN
                     WRITE( FMAT, '(A,I2,A)' )
     &                      '( ''+'', T', 24+(J-1)*13, ', 1P,E13.4 )'
                     WRITE( *, FMAT ) FNORM * PMOM( M, J )
                  END IF

   20          CONTINUE
   30       CONTINUE

         END IF

      END IF


      RETURN

      END

      SUBROUTINE SMALL1( XX, NUMANG, XMU, QEXT, QSCA, GQSC, SFORW,
     &                   SBACK, S1, S2, TFORW, TBACK, A, B )

c       Small-particle limit of Mie quantities in totally reflecting
c       limit ( Mie series truncated after 2 terms )

c        A,B       First two Mie coefficients, with numerator and
c                  denominator expanded in powers of XX ( a factor
c                  of XX**3 is missing but is restored before return
c                  to calling program )  ( Ref. 2, p. 1508 )

      IMPLICIT  NONE

c     .. Parameters ..

      DOUBLE PRECISION
     &   FIVNIN,FIVTHR,HALF,ONE,THRHLF,TWOTHR,TWO,ZERO
      PARAMETER ( 
     &   ZERO = 0.D0,
     &   HALF = 0.5D0,
     &   FIVNIN = 5.D0/9.D0,
     &   TWOTHR = 2.D0/3.D0,
     &   ONE = 1.D0,
     &   FIVTHR = 5.D0/3.D0,
     &   TWO = 2.D0,
     &   THRHLF = 1.5D0 )
c     ..
c     .. Scalar Arguments ..

      INTEGER
     &   NUMANG
      DOUBLE PRECISION
     &   GQSC, QEXT, QSCA, XX
      DOUBLE COMPLEX
     &   SBACK, SFORW
c     ..
c     .. Array Arguments ..

      DOUBLE PRECISION
     &   XMU( * )
      DOUBLE COMPLEX
     &   A( * ), B( * ), S1( * ), S2( * ), TBACK( * ), TFORW( * )
c     ..
c     .. Local Scalars ..

      INTEGER J
      DOUBLE PRECISION RTMP
      DOUBLE COMPLEX CTMP
c     ..
c     .. Intrinsic Functions ..

c 01.11.05 (BTD) attempt to make portable:
c      INTRINSIC DIMAG, CONJG, DCMPLX, DBLE, REALPART
      INTRINSIC DIMAG, DCONJG, DCMPLX, DBLE
c     ..
c     .. Statement Functions ..

      DOUBLE PRECISION SQ
      DOUBLE COMPLEX CJG
c     ..
c     .. Statement Function definitions ..

c 01.11.05 (BTD) attempt to make portable
c g77 is unhappy with following statement
c      SQ( CTMP ) = DBLE( CTMP )**2 + DIMAG( CTMP )**2
c solaris f77 is unhappy with following statement
c      SQ( CTMP ) = REALPART( CTMP )**2 + DIMAG( CTMP )**2
c try following:

      CJG(CTMP)=DCONJG(CTMP)
      SQ(CTMP)=ABS(CTMP)**2
c     ..

c                                                       ** Eq. R/A.5
c 01.11.05 (BTD) replace COMPLEX by DCMPLX
c      A( 1 ) = COMPLEX( ZERO, TWOTHR*( ONE - 0.2D0*XX**2 ) ) /
c     &         COMPLEX( ONE - HALF*XX**2, TWOTHR*XX**3 )
c                                                       ** Eq. R/A.6
c      B( 1 ) = COMPLEX( ZERO, - ( ONE - 0.1D0*XX**2 ) / 3.D0) /
c     &         COMPLEX( ONE + HALF*XX**2, - XX**3 / 3.D0)
c                                                       ** Eq. R/A.7,8
c      A( 2 ) = COMPLEX( ZERO,   XX**2 / 30.D0)
c      B( 2 ) = COMPLEX( ZERO, - XX**2 / 45.D0)

      A( 1 ) = DCMPLX( ZERO, TWOTHR*( ONE - 0.2D0*XX**2 ) ) /
     &         DCMPLX( ONE - HALF*XX**2, TWOTHR*XX**3 )
c                                                       ** Eq. R/A.6
      B( 1 ) = DCMPLX( ZERO, - ( ONE - 0.1D0*XX**2 ) / 3.D0) /
     &         DCMPLX( ONE + HALF*XX**2, - XX**3 / 3.D0)
c                                                       ** Eq. R/A.7,8
      A( 2 ) = DCMPLX( ZERO,   XX**2 / 30.D0)
      B( 2 ) = DCMPLX( ZERO, - XX**2 / 45.D0)
c                                                       ** Eq. R/A.9
      QSCA = 6.D0* XX**4 *( SQ( A(1) ) + SQ( B(1) ) +
     &           FIVTHR*( SQ( A(2) ) + SQ( B(2) ) ) )
      QEXT = QSCA
c                                                       ** Eq. R/A.10
      GQSC = 6.D0* XX**4 *DBLE( A(1)*CJG( A(2) + B(1) ) +
     &         ( B(1) + FIVNIN*A(2) )*CJG( B(2) ) )

      RTMP   = THRHLF * XX**3
      SFORW  = RTMP*( A(1) + B(1) + FIVTHR*( A(2) + B(2) ) )
      SBACK  = RTMP*( A(1) - B(1) - FIVTHR*( A(2) - B(2) ) )
      TFORW( 1 ) = RTMP*( B(1) + FIVTHR*( TWO*B(2) - A(2) ) )
      TFORW( 2 ) = RTMP*( A(1) + FIVTHR*( TWO*A(2) - B(2) ) )
      TBACK( 1 ) = RTMP*( B(1) - FIVTHR*( TWO*B(2) + A(2) ) )
      TBACK( 2 ) = RTMP*( A(1) - FIVTHR*( TWO*A(2) + B(2) ) )

      DO 10 J = 1, NUMANG
c                                                    ** Eq. R/A.11,12

         S1( J ) = RTMP*( A(1) + B(1)*XMU( J ) +
     &                    FIVTHR*( A(2)*XMU( J ) + 
     &                             B(2)*( TWO*XMU( J )**2 - ONE) ) )
         S2( J ) = RTMP*( B(1) + A(1)*XMU( J ) +
     &                    FIVTHR*( B(2)*XMU( J ) + 
     &                             A(2)*( TWO*XMU( J )**2 - ONE) ) )
   10 CONTINUE

c                                     ** Recover actual Mie coefficients
      A( 1 ) = XX**3 * A(1)
      A( 2 ) = XX**3 * A(2)
      B( 1 ) = XX**3 * B(1)
      B( 2 ) = XX**3 * B(2)

      RETURN
      END

      SUBROUTINE SMALL2( XX, CIOR, CALCQE, NUMANG, XMU, QEXT, QSCA,
     &                   GQSC, SFORW, SBACK, S1, S2, TFORW, TBACK,
     &                   A, B )

c       Small-particle limit of Mie quantities for general refractive
c       index ( Mie series truncated after 2 terms )

c        A,B       First two Mie coefficients, with numerator and
c                  denominator expanded in powers of XX ( a factor
c                  of XX**3 is missing but is restored before return
c                  to calling program )

c        CIORSQ    Square of refractive index


      IMPLICIT  NONE

c     .. Parameters ..

      DOUBLE PRECISION
     &   FIVTHR,HALF,ONE,THRHLF,TWO,TWOTHR,ZERO
      PARAMETER (
     &   ZERO = 0.D0,
     &   HALF = 0.5D0,
     &   TWOTHR = 2.D0/3.D0,
     &   ONE = 1.D0,
     &   THRHLF = 1.5D0,
     &   FIVTHR = 5.D0/3.D0,
     &   TWO = 2.D0)
      DOUBLE COMPLEX
     &   CXZERO
      PARAMETER(
     &   CXZERO=(0.D0,0.D0))
c     ..
c     .. Scalar Arguments ..

      LOGICAL
     &   CALCQE
      INTEGER
     &   NUMANG
      DOUBLE PRECISION
     &   GQSC, QEXT, QSCA, XX
      DOUBLE COMPLEX
     &   CIOR, SBACK, SFORW
c     ..
c     .. Array Arguments ..

      DOUBLE PRECISION
     &   XMU( * )
      DOUBLE COMPLEX
     &   A( * ), B( * ), S1( * ), S2( * ), TBACK( * ), TFORW( * )
c     ..
c     .. Local Scalars ..

      INTEGER
     &   J
      DOUBLE PRECISION
     &   RTMP
      DOUBLE COMPLEX
     &   CIORSQ, CTMP
c     ..
c     .. Intrinsic Functions ..

c 01.11.05 (BTD) attempt to make portable
c      INTRINSIC DIMAG, CONJG, DCMPLX, DBLE, REALPART
      INTRINSIC DIMAG, DCONJG, DCMPLX, DBLE
c     ..
c     .. Statement Functions ..

      DOUBLE PRECISION
     &   SQ
      DOUBLE COMPLEX
     &   CJG
c     ..
c     .. Statement Function definitions ..

c 01.11.05 (BTD) attempt to make portable
c g77 is unhappy with following statement
c      SQ( CTMP ) = DBLE( CTMP )**2 + DIMAG( CTMP )**2
c solaris f77 is unhappy with following statement
c      SQ( CTMP ) = REALPART( CTMP )**2 + DIMAG( CTMP )**2
c try following:

      CJG(CTMP)=DCONJG(CTMP)
      SQ(CTMP)=ABS(CTMP)**2
c     ..


      CIORSQ = CIOR**2
c 01.11.05 (BTD) replace COMPLEX by DCMPLX
c      CTMP   = COMPLEX( ZERO, TWOTHR )*( CIORSQ - ONE )
      CTMP   = DCMPLX( ZERO, TWOTHR )*( CIORSQ - ONE )

c                                           ** Eq. R42a
      A( 1 ) = CTMP*( ONE- 0.1D0*XX**2 +
     &         ( CIORSQ / 350.D0 + ONE/280.D0)*XX**4 ) /
     &         ( CIORSQ + TWO+ ( ONE- 0.7D0*CIORSQ )*XX**2 -
     &         ( CIORSQ**2 / 175.D0- 0.275D0*CIORSQ + 0.25D0 )*XX**4 +
     &         XX**3 * CTMP * ( ONE- 0.1D0*XX**2 ) )

c                                           ** Eq. R42b
      B( 1 ) = ( XX**2 / 30.D0 )*CTMP*( ONE+
     &         ( CIORSQ / 35.D0 - ONE/ 14.D0)*XX**2 ) /
     &         ( ONE- ( CIORSQ / 15.D0 - ONE/6.D0)*XX**2 )

c                                           ** Eq. R42c

      A( 2 ) = ( 0.1D0*XX**2 )*CTMP*( ONE- XX**2 / 14.D0 ) /
     &         ( TWO*CIORSQ + 3.D0- ( CIORSQ / 7.D0- HALF ) * XX**2 )

c                                           ** Eq. R40a

      QSCA = (6.D0*XX**4) * ( SQ( A(1) ) + SQ( B(1) ) +
     &                     FIVTHR * SQ( A(2) ) )

c                                           ** Eq. R40b
      QEXT = QSCA
      IF( CALCQE ) QEXT = 6.D0*XX * DBLE( A(1) + B(1) + FIVTHR*A(2) )

c                                           ** Eq. R40c

      GQSC = (6.D0*XX**4) * DBLE( A(1)*CJG( A(2) + B(1) ) )

      RTMP   = THRHLF * XX**3
      SFORW  = RTMP*( A(1) + B(1) + FIVTHR*A(2) )
      SBACK  = RTMP*( A(1) - B(1) - FIVTHR*A(2) )
      TFORW( 1 ) = RTMP*( B(1) - FIVTHR*A(2) )
      TFORW( 2 ) = RTMP*( A(1) + TWO*FIVTHR*A(2) )
      TBACK( 1 ) = TFORW(1)
      TBACK( 2 ) = RTMP*( A(1) - TWO*FIVTHR*A(2) )


      DO 10 J = 1, NUMANG
c                                      ** Eq. R40d,e

         S1( J ) = RTMP*( A(1) + ( B(1) + FIVTHR*A(2) )*XMU( J ) )
         S2( J ) = RTMP*( B(1) + A(1)*XMU( J ) +
     &                    FIVTHR*A(2)*( TWO*XMU( J )**2 - ONE) )
   10 CONTINUE

c                                     ** Recover actual Mie coefficients
      A( 1 ) = XX**3 * A(1)
      A( 2 ) = XX**3 * A(2)
      B( 1 ) = XX**3 * B(1)
      B( 2 ) = CXZERO

      RETURN
      END

      SUBROUTINE TESTMI( COMPAR, XX, CREFIN, MIMCUT, PERFCT, ANYANG,
     &                   NMOM, IPOLZN, NUMANG, XMU, QEXT, QSCA, GQSC,
     &                   SFORW, SBACK, S1, S2, TFORW, TBACK, PMOM,
     &                   MOMDIM )

c         Set up to run test case when  COMPAR = False;  when  = True,
c         compare Mie code test case results with correct answers
c         and abort if even one result is inaccurate.

c         The test case is :  Mie size parameter = 10
c                             refractive index   = 1.5 - 0.1 i
c                             scattering angle = 140 degrees
c                             1 Sekera moment

c         Results for this case may be found among the test cases
c         at the end of reference (1).

c         *** NOTE *** When running on some computers, esp. in single
c         precision, the Accur criterion below may have to be relaxed.
c         However, if Accur must be set larger than 10**-3 for some
c         size parameters, your computer is probably not accurate
c         enough to do Mie computations for those size parameters.

c     Routines called :  ERRMSG, MIPRNT, TSTBAD

      IMPLICIT  NONE

c     .. Scalar Arguments ..

      LOGICAL
     &   ANYANG, COMPAR, PERFCT
      INTEGER
     &   IPOLZN, MOMDIM, NMOM, NUMANG
      DOUBLE PRECISION
     &   GQSC, MIMCUT, QEXT, QSCA, XX
      DOUBLE COMPLEX
     &   CREFIN, SBACK, SFORW

c     .. Array Arguments ..

      DOUBLE PRECISION
     &   PMOM( 0:MOMDIM, * ), XMU( * )
      DOUBLE COMPLEX
     &   S1( * ), S2( * ), TBACK( * ), TFORW( * )

c     .. Local Scalars ..

      LOGICAL
     &   ANYSAV, OK, PERSAV
      INTEGER
     &   IPOSAV, M, N, NMOSAV, NUMSAV
      DOUBLE PRECISION
     &   ACCUR, CALC, EXACT, MIMSAV, TESTGQ, TESTQE, TESTQS,
     &   XMUSAV, XXSAV
      DOUBLE COMPLEX
     &   CRESAV, TESTS1, TESTS2, TESTSB, TESTSF

c     .. Local Arrays ..

      LOGICAL   CALCMO( 4 ), PRNT( 2 )
      DOUBLE PRECISION
     &   TESTPM( 0:1 )
      DOUBLE COMPLEX
     &   TESTTB( 2 ), TESTTF( 2 )

c     .. variables to be SAVEd ..

      SAVE      XXSAV, CRESAV, MIMSAV, PERSAV, ANYSAV, NMOSAV, IPOSAV,
     &          NUMSAV, XMUSAV

c     .. External Functions ..

      LOGICAL   TSTBAD
      EXTERNAL  TSTBAD

c     .. External Subroutines ..

      EXTERNAL  ERRMSG, MIPRNT

c     .. Intrinsic Functions ..

      INTRINSIC ABS, DIMAG, DBLE

c     .. Statement Functions ..

      LOGICAL WRONG
      WRONG( CALC, EXACT ) = ABS( ( CALC - EXACT ) / EXACT ).GT.ACCUR

      DATA      TESTQE / 2.459791 /,
     &          TESTQS / 1.235144 /,
     &          TESTGQ / 1.139235 /,
     &          TESTSF / ( 61.49476, -3.177994 ) /,
     &          TESTSB / ( 1.493434, 0.2963657 ) /,
     &          TESTS1 / ( -0.1548380, -1.128972 ) /,
     &          TESTS2 / ( 0.05669755, 0.5425681 ) /,
     &          TESTTF / ( 12.95238, -136.6436 ),
     &                   ( 48.54238, 133.4656 ) /,
     &          TESTTB / ( 41.88414, -15.57833 ),
     &                   ( 43.37758, -15.28196 ) /,
     &          TESTPM / 227.1975, 183.6898 /

      DATA      ACCUR / 1.D-4 /
c     ..


      IF( .NOT.COMPAR ) THEN
c                                   ** Save certain user input values
         XXSAV  = XX
         CRESAV = CREFIN
         MIMSAV = MIMCUT
         PERSAV = PERFCT
         ANYSAV = ANYANG
         NMOSAV = NMOM
         IPOSAV = IPOLZN
         NUMSAV = NUMANG
         XMUSAV = XMU( 1 )
c                                   ** Reset input values for test case
         XX     = 10.D0
         CREFIN = ( 1.5D0, -0.1D0 )
         MIMCUT = 0.D0
         PERFCT = .FALSE.
         ANYANG = .TRUE.
         NMOM   = 1
         IPOLZN = -1
         NUMANG = 1
         XMU( 1 ) = -0.7660444D0

      ELSE
c                                    ** Compare test case results with
c                                    ** correct answers and abort if bad
         OK = .TRUE.

         IF( WRONG( QEXT,TESTQE ) )
     &       OK = TSTBAD( 'QEXT', ABS( ( QEXT - TESTQE ) / TESTQE ) )

         IF( WRONG( QSCA,TESTQS ) )
     &       OK = TSTBAD( 'QSCA', ABS( ( QSCA - TESTQS ) / TESTQS ) )

         IF( WRONG( GQSC,TESTGQ ) )
     &       OK = TSTBAD( 'GQSC', ABS( ( GQSC - TESTGQ ) / TESTGQ ) )

         IF( WRONG( DBLE( SFORW ),DBLE( TESTSF ) ) .OR.
     &       WRONG( DIMAG( SFORW ),DIMAG( TESTSF ) ) )
     &       OK = TSTBAD( 'SFORW', ABS( ( SFORW - TESTSF ) / TESTSF ) )

         IF( WRONG( DBLE( SBACK ),DBLE( TESTSB ) ) .OR.
     &       WRONG( DIMAG( SBACK ),DIMAG( TESTSB ) ) )
     &       OK = TSTBAD( 'SBACK', ABS( ( SBACK - TESTSB ) / TESTSB ) )

         IF( WRONG( DBLE( S1(1) ),DBLE( TESTS1 ) ) .OR.
     &       WRONG( DIMAG( S1(1) ),DIMAG( TESTS1 ) ) )
     &       OK = TSTBAD( 'S1', ABS( ( S1(1) - TESTS1 ) / TESTS1 ) )

         IF( WRONG( DBLE( S2(1) ),DBLE( TESTS2 ) ) .OR.
     &       WRONG( DIMAG( S2(1) ),DIMAG( TESTS2 ) ) )
     &       OK = TSTBAD( 'S2', ABS( ( S2(1) - TESTS2 ) / TESTS2 ) )


         DO 10  N = 1, 2

            IF( WRONG( DBLE( TFORW(N) ),DBLE( TESTTF(N) ) ) .OR.
     &          WRONG( DIMAG( TFORW(N) ),
     &          DIMAG( TESTTF(N) ) ) ) OK = TSTBAD( 'TFORW',
     &          ABS( ( TFORW(N) - TESTTF(N) ) / TESTTF(N) ) )

            IF( WRONG( DBLE( TBACK(N) ),DBLE( TESTTB(N) ) ) .OR.
     &          WRONG( DIMAG( TBACK(N) ),
     &          DIMAG( TESTTB(N) ) ) ) OK = TSTBAD( 'TBACK',
     &          ABS( ( TBACK(N) - TESTTB(N) ) / TESTTB(N) ) )

   10    CONTINUE

         DO 20 M = 0, 1

            IF ( WRONG( PMOM(M,1), TESTPM(M) ) )
     &           OK =  TSTBAD( 'PMOM', ABS( (PMOM(M,1)-TESTPM(M)) /
     &                                      TESTPM(M) ) )

   20    CONTINUE

         IF( .NOT.OK ) THEN

            PRNT( 1 ) = .TRUE.
            PRNT( 2 ) = .TRUE.
            CALCMO( 1 ) = .TRUE.
            CALCMO( 2 ) = .FALSE.
            CALCMO( 3 ) = .FALSE.
            CALCMO( 4 ) = .FALSE.

            CALL MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT,
     &                   QSCA, GQSC, NMOM, IPOLZN, MOMDIM, CALCMO, PMOM,
     &                   SFORW, SBACK, TFORW, TBACK, S1, S2 )

            CALL ERRMSG( 'MIEV0 -- Self-test failed',.TRUE.)

         END IF
c                                       ** Restore user input values
         XX     = XXSAV
         CREFIN = CRESAV
         MIMCUT = MIMSAV
         PERFCT = PERSAV
         ANYANG = ANYSAV
         NMOM   = NMOSAV
         IPOLZN = IPOSAV
         NUMANG = NUMSAV
         XMU( 1 ) = XMUSAV

      END IF

      RETURN
      END
      SUBROUTINE MSPHSX(X,MU,EPS,NTHETA,THETA,
     &                  QEXT,QSCA,QABS,QBAC,S1,S2)
c
c **************************************************
c
c Subroutine msphsx computes an approximate values of the scattered amplitude!
c and scattering cross sections (efficiencies) for a magnetic sphere with
c small size parameters.
c
c     Merrill Milham        >>> version 1.2 <<<             DEC 1993
c
c     Inputs :
c              x = size parameter of the sphere (real*8)
c            eps = complex permittivity: epsr -i*epsi (complex*16)
c             mu = complex permeability: mur - i*mui (complex*16)
c         ntheta = number of scattering angles (integer)
c          theta = scattering angles in degrees (real*8)
c
c Outputs :
c             s1 = scattered amplitude
c             s2 = scattered amplitude
c           qext = extinction efficiency
c           qsca = scattering efficiency
c           qabs = absorption efficiency
c           qbac = backscatter efficiency
c
c References:
c
c M. Kerker, D.-S. Wang, and C. L. Giles, "Electromagnetic scattering
c by magnetic spheres," J. Opt. Soc. Am., 73, 765-767 (1983).
c
c J. A. Stratton, "Electromagnetic Theory," (McGraw-Hill, New York,
c 1941)
c
c ***********************************************
c
      REAL*8 X
      COMPLEX*16 MU,EPS
      INTEGER NTHETA     
      REAL*8 THETA(1)
*
      COMPLEX*16 S1(1),S2(1)
      REAL*8 QEXT,QSCA,QABS,QBAC
*
	REAL*8 XP,N,K,THETAL
      COMPLEX*16 EF,MF,MUI,EPSI,CX
*
      REAL*8 ZERO,ONE,TWO,THREE,FOUR,EIGHT,C0,CPI,RAD
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,FOUR=4.D0)
      PARAMETER (THREE=3.D0,EIGHT=8.D0,C0=EIGHT/THREE)
      PARAMETER (CPI=3.1415926535897932384D0,RAD=CPI/180.D0)
      COMPLEX*16 CDBLE0,CDBLEI
c
c page 34
c
      PARAMETER (CDBLE0=(0.D0,0.D0),CDBLEI=(0.D0,1.D0)) 
*
      INTEGER NANGL,NANGL2
      PARAMETER (NANGL=255,NANGL2=(NANGL+1)/2)
cthw      real*8 xmu(nangl)
      REAL*16 XMU(NANGL)
      INTEGER L,J,J2,J3
c
      IF(NTHETA.EQ.0) THEN
            S1(1)=CDBLE0
            S2(1)=CDBLE0
                        ELSE
*
                IF(NTHETA.LE.NANGL2) THEN
                    CONTINUE
                                      ELSE
       PRINT *,NTHETA,'scattering angles input: only',NANGL2,' allowed'
                    STOP
                ENDIF 
*
      DO 100 L=1,NTHETA
      THETAL=DABS(THETA(L))
      THETA(L)=THETAL 
*
               IF(THETAL.LE.90.D0) THEN
                CONTINUE
                                    ELSE
      PRINT*,'theta(',l,')=',THETAL,'scattering angles must be <=90 deg'
            STOP
            END IF 
*
      THETAL=RAD*THETAL
      XMU(L)=DCOS(THETAL)
100   CONTINUE
*
      END IF
c
      N=DBLE(EPS)
      K=DABS(DIMAG(EPS))
      EPSI=DCMPLX(N,-K) 
*
      N=DBLE(MU)
      K=DABS(DIMAG(MU))
      MUI=DCMPLX(N,-K)
c
      EF=(EPSI-ONE)/(EPSI+TWO)
      MF=(MUI-ONE)/(MUI+TWO) 
*
      XP=X*X*X
      CX=DCMPLX(ZERO,XP)
c
      IF(NTHETA.EQ.0) THEN
            CONTINUE
c
c page 35
c
            ELSE
      	J2=2*NTHETA
      	DO 500 J=1,NTHETA
      	J3=J2-J
      	S1(J)=CX*(EF+MF*XMU(J))
      	S2(J)=CX*(EF*XMU(J)+MF)
      	S1(J3)=CX*(EF-MF*XMU(J))
      	S2(J3)=CX*(MF-EF*XMU(J))
500   CONTINUE
      END IF
c
      XP=XP*X
      QEXT=-FOUR*X*DIMAG(EF+MF)
cthw	  qsca=c0*xp*(zabs(ef)**2+(mf)**2)
      QSCA=C0*XP*(ABS(EF)**2+ABS(MF)**2)
      
      IF(QEXT.LE.ZERO) QEXT=QSCA
      QABS=DDIM(QEXT,QSCA)
cthw      qbac=four*xp*(zabs(ef-mf))**2
      QBAC=FOUR*XP*(ABS(EF-MF))**2
c
      RETURN
c
      END
      SUBROUTINE PAH_CRS_SCT(NC,NH,D2H,ZG,WAVELENGTHCM,CROSS_SECTION)
      IMPLICIT NONE

c=============== pah_crs_sct_v15 ====================================
c arguments:

      INTEGER NC,NH,ZG
      DOUBLE PRECISION CROSS_SECTION,D2H,WAVELENGTHCM

c common:
c /FIP_COM/ is for communication between pah_crs_sct and qcomp
c                  (and from qcomp to jnu_u)

      DOUBLE PRECISION FIP
      COMMON/FIP_COM/FIP

c local variables:

c*** diagnostic
c      integer iabort
c***
      CHARACTER STATE*7
      DOUBLE PRECISION 
     &   C1,C2,C3,C4,CC_SHIFT,CDBEND_SHIFT,CD_REDUCED,CDSTRE_SHIFT,
     &   CRS_IP,CRS_SCT_FEATURE,CRS_SCT_FIRCONT,CRS_SCT_NIR,CRS_SCT_UV,
     &   CS33,CS_VISUAL,CX2XC,CX2XC33,
     &   DCRS,
     &   ENHANCE_033I,ENHANCE_033N,ENHANCE_052,ENHANCE_062,ENHANCE_077,
     &   ENHANCE_086,ENHANCE_113,ENHANCE_120,ENHANCE_126,ENHANCE_135,
     &   ENHANCE_142,ENHANCE_159,ENHANCE_164,ENHANCE_17,ENHANCE_189,
     &   F1,F2,F3,FIP0722,FIP2175,FIP_033,FIP_052,FIP_057,FIP_062,
     &   FIP_067,FIP_074,FIP_076,FIP_079,FIP_083,FIP_086,FIP_107,
     &   FIP_112,FIP_113,FIP_120,FIP_126,FIP_127,FIP_135,FIP_142,
     &   FIP_159,FIP_164,FIP_17,FIP_170,FIP_174,FIP_179,FIP_189,FIP_FIR,
     &   FRACT_D,FRACT_H,
     &   GAMMA,GAMMA0722,GAMMA2175,
     &   GAMMA_033,GAMMA_052,GAMMA_057,GAMMA_062,GAMMA_067,GAMMA_074,
     &   GAMMA_076,GAMMA_079,GAMMA_083,GAMMA_086,GAMMA_107,GAMMA_112,
     &   GAMMA_113,GAMMA_120,GAMMA_126,GAMMA_127,GAMMA_135,GAMMA_142,
     &   GAMMA_159,GAMMA_164,GAMMA_170,GAMMA_174,GAMMA_179,GAMMA_189,
     &   GAMMA_FIR,
     &   H2C,
     &   INT_CROSSCT0722,INT_CROSSCT2175,
     &   INT_CROSSCT_033,INT_CROSSCT_052,INT_CROSSCT_057,
     &   INT_CROSSCT_062,INT_CROSSCT_067,INT_CROSSCT_074,
     &   INT_CROSSCT_076,INT_CROSSCT_079,INT_CROSSCT_083,
     &   INT_CROSSCT_086,INT_CROSSCT_107,INT_CROSSCT_112,
     &   INT_CROSSCT_113,INT_CROSSCT_120,INT_CROSSCT_126,
     &   INT_CROSSCT_127,INT_CROSSCT_135,INT_CROSSCT_142,
     &   INT_CROSSCT_159,INT_CROSSCT_164,INT_CROSSCT_17,
     &   INT_CROSSCT_170,INT_CROSSCT_174,INT_CROSSCT_179,
     &   INT_CROSSCT_189,INT_CROSSCT_FIR,INT_CROSSCTD_033,
     &   INT_CROSSCTD_062,INT_CROSSCTD_074,INT_CROSSCTD_076,
     &   INT_CROSSCTD_079,INT_CROSSCTD_083,INT_CROSSCTD_086,
     &   INT_CROSSCTD_112,INT_CROSSCTD_113,INT_CROSSCTD_120,
     &   INT_CROSSCTD_126,INT_CROSSCTD_127,NRING,
     &   WAVEUM,WIDTH1,WIDTH2,
     &   WL0722,WL2175,
     &   WL_033,WL_052,WL_057,WL_062,WL_067,WL_074,WL_076,WL_079,WL_083,
     &   WL_086,WL_107,WL_112,WL_113,WL_120,WL_126,WL_127,WL_135,WL_142,
     &   WL_159,WL_164,WL_170,WL_174,WL_179,WL_189,WL_FIR,
     &   WLD_033,WLD_062,WLD_074,WLD_076,WLD_079,WLD_083,WLD_086,
     &   WLD_112,WLD_113,WLD_120,WLD_126,WLD_127,
     &   X,X2XC,XC,XPEAK1,XPEAK2

      DOUBLE PRECISION DRUDES
      EXTERNAL DRUDES
c ---------------------------------------------------------------------
c Subroutine PAH_CRS_SCT
c Purpose: calculates the UV/IR cross section of neutral/ionized PAHs.        -
c input:    
c       NC -- the number of carbon atoms (integer);  
c       NH -- the number of hydrogen atoms (integer; sum of H and D);  
c       D2H -- deuterium to total hydrogen ratio = D/(H+D) ; 
c       ZG -- grain charge (integer); 
c             currently, only consider zg=0 (neutral)
c             or zg!=0 (PAH cation or anion);
c       WAVELENGTHCM -- in cm!

c output:
c       CROSS_SECTION = absorption cross section per C atom 
c                       (cm^2/c-atom) for randomly-oriented PAH
c       FIP = fraction of absorption contributed by in-plane vibrations
c-----------------------------------------------------------------------
c history
c originally written by Aigen Li, Princeton University Observatory
c 1999.12.22 (AL)  no dehydrogenation effects were considered;
c 1999.12.28 (AL)  when modeling the interstellar IR emission, 
c                  the strength of 6.2, 7.7, 8.6 micron features need be 
c                  enhanced by a factor of 1, 2, or 3 (enhancement).
c 2000.03.08 (AL)  add deuterated PAH's emission features -- 
c                  upon deuteration, CH shifts to lower frequency by 1.2-1.3; 
c                  C-C by 1.1 (Hudgins et al. 94); 
c                  The C-H strength reduced by 1.75 (Bauschlicher et al. 
c                  1997). The C-H stretching mode shifts by 1.35; 
c                  C-H bending mode shifts 1.25; C-C mode no shifts
c 2005.10.31 (BTD) cosmetic changes in preparation for modifications
c 2005.10.31 (BTD) substantial rewrite
c                  have PAH emission features have been modified
c                  (but this is preliminary -- additional modifications
c                  as informed by SINGS spectra).
c                  * adjusted wavelength and gamma for 6.2um feature
c                    wave  6.2   -> 6.22
c                    gamma 0.032 -> 0.0284
c                  * resolved 7.7um feature into 3 components
c                    wave   gamma p/ptot
c                    7.417  0.126 0.26
c                    7.598  0.044 0.39
c                    7.850  0.053 0.35
c                  * resolved 8.6um feature into 2 components
c                    wave   gamma p/ptot
c                    8.329  0.052 0.58
c                    8.610  0.039 0.42
c                  * resolved 11.3um feature into 2 components
c                    11.23  0.010 0.18
c                    11.29  0.029 0.82
c                  * adjusted wavelength and gamma for 11.9um feature
c                    wave  11.9 -> 12.00
c                    gamma .025 -> 0.050
c                  * adjusted wavelength and gamma for 12.6um feature
c                    wave  12.7 -> 12.60
c                    gamma .024 -> 0.043
c                  * adjusted wavelength and gamma for 16.4um feature
c                    wave  16.4 -> 16.447
c                    gamma 0.01 -> 0.014
c                  * add 17um complex
c                     wave  gamma p/ptot
c                    17.038 0.065 0.50
c                    17.377 0.012 0.25
c                    17.873 0.016 0.25
c                    try INT_CROSSCT_17=10.D-20 cm
c                  * removed features at 18.3, 21.2, 23.1um
c                  * modify FIR continuum component
c                    sigma_int = 15.2d-20 cm -> 20.d-20 cm
c                    wave        26um        -> 20um
c                    gamma       0.69        -> 0.7
c 2005.11.28 (BTD) * modify params for 8.33 and 8.61um features
c                    (had incorrect strength ratio)
c                  * add features at 13.6 and 14.3um based on SINGS
c                    spectra
c                  * modify params for FIR continuum component
c                    wave        20um        -> 17um
c                    gamma       0.70        -> 1.0
c 2005.11.29 (BTD) * add NIR absorption term for pah ions following
c                    Mattioda, Allamandola, & Hudgins 2005, 
c                    ApJ 629, 1183
c 2005.12.05 (BTD) * modify params for 17.04,17.38,17.87um features
c                    (had incorrect strength ratio)
c                    new values: 
c                     wave  gamma p/ptot
c                    17.038 0.065 0.90
c                    17.377 0.012 0.06
c                    17.873 0.016 0.04
c 2005.12.12 (BTD) * Adjusted Mattioda et al fit parameters for NIR
c                    absorption
c 2006.04.15 (BTD) * Enhanced strength of 16.45,17.04,17.38,17.87 
c                    by factor ENHANCE_17=2.
c                  * Shifted PAH broad drude component from
c                    w0=17um,gamma=0.7,sigint=5e-19cm
c                    to
c                    w0=15um,
c                    gamma=0.8,
c                    sigint=7e-19*ENHANCE_17=1.4e-18cm 
c 2006.05.09 (BTD) * v4: Enhance strength of 3.3um feature by factor 
c                    ENHANCE_033N=2.0 for neutrals,
c                    ENHANCE_033I=3.0 for anions and cations
c 2006.05.18 (BTD) * v5: Reduce strength of neutral 3.3um feature to
c                    ENHANCE_033N=1.5 for neutrals
c                    ENHANCE_033I=2.0 for anions and cations
c 2006.06.07 (BTD) * v6: Add features at 5.70, 6.69, 15.90um
c                    Change 14.30 -> 14.19um
c                    change "143" -> "142"
c                    Change interpretation of input D2H
c                    previous: D2H = ^2H/^1H  : 0 .le. D2H .le. infty
c                    now: D2H = ^2H/(^1H+^2H) : 0 .le. D2H .le. 1
c                    New View of Deuteration:
c                    Because we are now using a detailed representation
c                    of observed 5-20um spectra based on astronomical
c                    objects that ARE partially deuterated, it does not
c                    make sense to try to add deuterated features
c                    *except* for the 4.40um C-D stretch
c                    Have therefore disabled addition of deuterated
c                    emission features except for 4.40um C-D mode
c 2006.06.24 (BTD) * v7: adjust features in view of theoretical data
c                    neutral:
c                       raise 3.3um strength by factor 4/3
c                             ENHANCE_033N=1.5 -> ENHANCE_033N=2.
c                    ion and neutral:
c                       lower 6.2,7.7,8.6,11.3 by factor 0.5
c                             ENHANCE_062=3. -> ENHANCE_062=1.5
c                             ENHANCE_077=2. -> ENHANCE_077=1.
c                             ENHANCE_086=2. -> ENHANCE_086=1.
c                             introduce         ENHANCE_113=0.5
c                                               ENHANCE_120=0.5
c                             ENHANCE_126=1.25->ENHANCE_126=0.625
c                             introduce         ENHANCE_136=0.5
c                                               ENHANCE_142=0.5
c                                               ENHANCE_159=0.5
c                             ENHANCE_17=2.  -> ENHANCE_17=1.
c 2006.06.25 (BTD) * v7: further adjustments
c                      ENHANCE_120=0.5 -> 1.0
c                      change INT_CROSSCT_067 to 0.2*INT_CROSSCT_062
c                      ENHANCE_170=1.0 -> 0.5
c                      add 5.25um feature with
c                         WL_052=5.25
c                         GAMMA_052=0.030
c                         INT_CROSSCT_052=20.D-20*ENHANCE_052 for ion
c                         INT_CROSSCT_052=2.5D-20*ENHANCE_052 for neutral
c                         ENHANCE_052=1 for first estimate
c 2006.09.07 (BTD) * v8
c                    5.27um feature: change 5.25um -> 5.27um
c                                    change gamma=0.030 -> 0.034
c                                    [Smith etal 2006 Table 4]
c                    5.70um feature: reduce gamma=0.040 -> 0.035
c                                    [Smith etal 2006 Table 4]
c                    6.22um feature: increase gamma=0.284 -> 0.030
c                                    [Smith etal 2006 Table 4]
c                    10.68um feature: add new feature using wavelength
c                                    and gamma from Smith etal 2006
c                                    and basing cross section on ratio of
c                                    10.68um power to 11.23um power
c                                    in SINGS sample
c                    11.23um feature: change gamma from 0.010 to 0.012
c                    11.33um feature: change wave from 11.30 to 11.33
c                                     change gamma from 0.029 to 0.032
c                    11.99um feature: change gamma from 0.050 to 0.045
c                    12.62um feature: change wave from 12.61 to 12.62
c                                     change gamma from .0435 to .042
c                    12.69um feature: add new feature using wavelength
c                                     and gamma from Smith etal 2006
c                                     and basing cross section on ratio
c                                     0.05=power in 12.69um feature/
c                                          power in 12.62um feature
c                    13.48um feature: change wave from 13.60 to 13.48
c                                     change gamma from 0.020 to 0.040
c                                     base cross section on ratio
c                                     0.30 = power in 13.48um feature/
c                                            power in 12.62um feature
c 06.10.24 (BTD) fixed numerical errors in values of
c                INT_CROSSCT_126 : 62.1 -> 31.0
c                                  69.6 -> 34.8
c                INT_CROSSCTD_126: 69.6 -> 34.8
c                INT_CROSSCT_142 : 1.0  -> 0.50
c                INT_CROSSCT_FIR : 100. -> 50.
c 06.10.28 (BTD) v9: 30% reduction in strength of 15-19um features:
c                ENHANCE_142     : 1.0 -> 0.9
c                ENHANCE_159       1.0 -> 0.8
c                ENHANCE_17        1.0 -> 0.7
c                ENHANCE_189       1.0 -> 0.7
c 06.10.31 (BTD) v9: redefine some variables to get round numbers
c                add new variable ENHANCE_164 = 1
c                set INT_CROSSCT_142 = 0.45D-20 (with ENHANCE_142=1)
c                set INT_CROSSCT_164 = 0.50D-20 (with ENHANCE_164=1)
c                set INT_CROSSCT_159 = 0.04D-20 (with ENHANCE_159=1)
c                set INT_CROSSCT_17  = 2.40D-20 (with ENHANCE_17=1)
c                set INT_CROSSCT_189 = 0.10D-20 (with ENHANCE_189=1)
c 06.11.01 (BTD) v9: corrected typos
c                WL113=11.30 -> 11.33
c                removed spurious (leftover) factor of 1/3 in 
c                        INT_CROSSCT_112 and INT_CROSSCT_113
c 07.12.10 (BTD) v9: changed ZG from DOUBLE PRECISION to INTEGER
c 08.01.04 (BTD) v10:
c                increase INT_CROSSCT_17  from 2.4D-20 to 3.5D-20
c                decrease INT_CROSSCT_FIR from 50.D-20 to 30.D-20
c 08.01.17 (BTD) v11:
c                increase INT_CROSSCT_17 from 3.5D-20 to 4.5D-20
c                decrease INT_CROSSCT_FIR from 30.D-2 to 25.D-20
c 08.03.17 (BTD) v12:
c                * increase INT_CROSSCT_142 from 0.45D-20 to 0.60D-20
c                * decrease INT_CROSSCT_17 from 4.5D-20 to 3.8D-20
c 08.03.27 (BTD) v12:
c                * discovered that we were inadvertently failing to
c                  include 12.69um feature ("127"): now corrected
c                * added 12.69um feature to deuterated calculation
c 08.03.31 (BTD) v13
c                * enhanced 8.3 and 8.6um features by 10%
c                  ENHANCE_086 = 1.0 -> 1.1
c                * weakend 17.um feature by 10% 
c                  ENHANCE_170 = 1.0 -> 0.9
c 08.04.08 (BTD) v14
c                * enhanced 8.3 and 8.6um features by 15%
c                  ENHANCE_086 = 1.1 -> 1.25
c                * weakened 17 um feature by 5%
c                  ENHANCE_170 = 0.9 -> 0.85
c 08.07.28 (BTD) v15
c                Added code to compute in-plane and out-of-plane fraction
c                for contribution to cross section at selected wavelength.
c                To accomplish this, introduce new parameters FIP_033, etc
c                identifing in-plane "fraction" for each mode.
c                For modes longward of 15.9um, we take FIP_xxx=(2./3.)
c                because vibrational character of these modes is not known.
c end history
c-----------------------------------------------------------------------

c*** diagnostic
c      iabort=0
c      if(wavelengthcm.gt.2.1e-5.and.wavelengthcm.lt.2.2e-5)then
c         write(0,*)'entered pah_cr_sct_v9 with NC,NH=',NC,NH,
c     &             ' WAVELENGTHCM=',WAVELENGTHCM
c      endif
c***

      IF(ZG.EQ.0)THEN
      STATE='neutral'	
      ELSEIF(ZG.NE.0)THEN
      STATE='ionized'
      ENDIF

      WAVEUM=WAVELENGTHCM*1.D4

c wavelengthcm: cm! wavelength: micron.
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c We first work on the UV (FUV-UV-visible) part.

      IF(NC.LE.40)NRING=0.3*NC
      IF(NC.GT.40)NRING=0.4*NC	

c nring: number of rings in a PAH with a number nc of carbon atoms. 

c C1,C2: linear part cros.sct=c1+c2/lambda;
c C4: FUV rise;      

      C1=1.8687D-18
      C2=0.1905D-18
      C4=0.7743D-18

      X=1.0/WAVEUM

c x: 1/micron

c WL0722    = wavelength adopted for sigma-sigma* resonance (um)
c GAMMA0722 = adopted gamma for sigma-sigma* resonance
c INT_CROSSCT0722 = integrated cross section (cm/C)

      WL0722=0.0722D0
      GAMMA0722=0.195D0
      INT_CROSSCT0722=7.97D-13

c 080723 BTD v15: sigma-sigma* resonance at 722nm=13.85um-1 is associated 
c                 with E perp c dielectric function.  Hence FIP0722=1.
      FIP0722=1.D0

c WL2175    = wavelength of "2175A feature" (um) 
c GAMMA2175 = adopted gamma for "2175A feature" (dimensionless)
c           = (.0473um FWHM)/(.2175um)
c INT_CROSSCT2175 = integrated cross section (cm/C)

      WL2175=0.2175D0
      GAMMA2175=0.2170D0
      INT_CROSSCT2175=1.23D-13

c 080723 BTD v15: 2175A feature is assumed to be due to pi-pi* transition 
c in aromatic C
c pi-pi* transition is associated with E perp c.  Hence FIP2175=1.
      FIP2175=1.D0

      CRS_IP=0.

      IF(X.GT.15.D0.AND.X.LE.17.25D0)THEN
         DCRS=126.D-18-6.4943D-18*X
         CRS_SCT_UV=DCRS

c 080723 BTD v15: at high frequencies, expect to have absorption from
c                 E parallel to c as well as E perp c.  Assume FIP=2/3
         CRS_IP=(2./3.)*DCRS

      ELSEIF(X.GT.10.D0.AND.X.LE.15.D0)THEN
         DCRS=-3.D-18+1.35D-18*X
         CRS_SCT_UV=DCRS
c 080723 BTD v15: take "continuum" to be 2/3 in-plane
         CRS_IP=(2./3.)*DCRS

         DCRS=INT_CROSSCT0722*DRUDES(WL0722,GAMMA0722,WAVEUM)
         CRS_SCT_UV=CRS_SCT_UV+DCRS
         CRS_IP=CRS_IP+FIP0722*DCRS

      ELSEIF(X.GT.7.7D0.AND.X.LE.10.D0)THEN
         DCRS=66.302D-18-24.367D-18*X+2.9501D-18*X*X-0.10569D-18*X**3
         CRS_SCT_UV=DCRS
c 080723 BTD v15: take "continuum" to be 2/3 in-plane
         CRS_IP=(2./3.)*DCRS

      ELSEIF(X.GT.5.9D0.AND.X.LE.7.7D0)THEN
         DCRS=C1+C2*X+C4*(0.5392*(X-5.9)**2+0.05644*(X-5.9)**3)+
     &        INT_CROSSCT2175*DRUDES(WL2175,GAMMA2175,WAVEUM)
         CRS_SCT_UV=DCRS
         CRS_IP=FIP2175*DCRS

      ELSEIF(X.GE.3.3D0.AND.X.LE.5.9D0)THEN
         DCRS=C1+C2*X+INT_CROSSCT2175*DRUDES(WL2175,GAMMA2175,WAVEUM)
         CRS_SCT_UV=DCRS
         CRS_IP=FIP2175*DCRS

      ELSEIF(X.LT.3.3D0)THEN

c PAH absorption edge differs from cation/anion to neutrals.

         IF(STATE.EQ.'ionized')THEN
            XC=1.D-4*(22824.*NRING**(-0.5)+8892.)
         ELSEIF(STATE.EQ.'neutral')THEN
            XC=1.D-4*(38040.*NRING**(-0.5)+10520.)
         ENDIF
         X2XC=X/XC

c CX2XC = cutoff function of Desert et al (1990):
c cabs is assumed to vary as 10.**(-3.431/X)*CX2XC

         CX2XC=ATAN(1000.*((X2XC-1.0)**3.0)/X2XC)/3.1416+0.5
         CX2XC33=ATAN(1000.*((3.3/XC-1.0)**3.0)/(3.3/XC))/3.1416+0.5

c normalize the near-UV/visible absorption to smoothly join 
c that of the Drude profile+linear continuum at 1/wave = 3.3um-1
c 2.6530 = .5392*(3.3-5.9)^2+.05644*(3.3-5.9)**3
c 3.431/3.3 = 1.03970

         CS33=C1+C2*3.3+
     &        INT_CROSSCT2175*DRUDES(WL2175,GAMMA2175,0.30303D0)
         DCRS=CS33*10.**(1.03970D0-3.431/X)*CX2XC/CX2XC33
         CRS_SCT_UV=DCRS

c 080723 BTD assume that near-UV/visible absorption is entirely in-plane
         CRS_IP=CRS_IP+DCRS

      ENDIF

c add NIR absorption for PAH ions as recommended by Mattioda,
c Allamandola, and Hudgins 2005
c NB: the drude term with lambda_3=1.915 is indeed negative
c     the factor exp(-(0.1/WAVEUM)**2) is introduced to cutoff
c     this absorption cross section at high energies.  This factor
c     is effectively unity for lambda > 1um where this NIR contribution
c     makes a difference.  Note that the integrated cross sections
c     given by Mattioda et al are in units of cm^2 um-1 C-1.
c     Here we need to convert to cm C-1: cm^2 um-1 C-1 = 1e4 cm C-1
c 2005.12.12 (BTD) A. Mattioda found error in fits and provided new
c                  values of Drude parameter
c                  lambda   gamma   S
c                  1.05um   0.055  2.0e-20 cm^2/um = 2.0e-16 cm C-1
c                  1.26     0.11   7.8e-21         = 7.8e-17 cm C-1
c                  1.905    0.09   1.465e-22       = 1.465e-18 cm C-1
c
      IF(STATE.EQ.'ionized')THEN
         CRS_SCT_NIR=3.5D-19*10.**(-1.45D0*WAVEUM)*
     &               EXP(-(0.1D0/WAVEUM)**2)+
     &               2.000D-16*DRUDES(1.050D0,0.055D0,WAVEUM)+
     &               7.800D-17*DRUDES(1.260D0,0.11D0,WAVEUM)-
     &               1.465D-18*DRUDES(1.905D0,0.09D0,WAVEUM)

      ELSE
         CRS_SCT_NIR=0.
      ENDIF

c 080723 BTD assume that NIR absorption is entirely in-plane:
      CRS_IP=CRS_IP+CRS_SCT_NIR

c      if(iabort.gt.0)then
c         cross_section=crs_sct_uv
c         return
c      endif

c crs_sct_uv: cm^2/c-atom.


c------------------- Calculation of IR/FIR absorption ---------------------
c This module calculates absorption due to discrete PAH vibrational modes
c from 3.3um and 18um
c plus an assumed continuum represented by a strongly-damped oscillator
c with central wavelength 20um and gamma=0.7
c
c Also includes possible deuteration of PAHs with resulting
c (1) shifts in C-H vibrational modes to C-D modes at longer wavelengths
c     fractional shift may be different for stretching and bending
c (2) reductions in strength of C-H -> C-D vibrational modes
c     reduction in strength allowed to be different for stretching and bending
c (3) allow for possible shift in C-C modes due to deuteration, but
c     at present shift is set to zero
c
c Relevant references:
c Pauzat et al. 1997, A&A, 319, 318: the effects of dehydrogenation
c on CH vibration features (3.3, 11.3, 11.9, 12.7);
c Allain et al. 1996, A&A, 305, 616: nc>30 fully hydrogenated; 
c deH: degree of dehydrogenation;
c
c-----------------------------------------------------------------------

c 06.06.07 (BTD) change CDSTRE_SHIFT from 1.35 to 1.330
c                (observed C-D stretching modes are
c                 at 4.356,4.389,4.397um in phenanthrene C14D10, 
c                    4.366,4.371,4.378,4.399,4.406,4.419 in pyrene C16D10
c                    4.395 in chyrsene C18D12
c                    mean of above 10 wavelengths = 4.389um
c                    4.389/3.300=1.330
c                
      CDSTRE_SHIFT=1.33D0
      CDBEND_SHIFT=1.25D0

      CC_SHIFT=1.00D0
      CD_REDUCED=1.D0/1.75D0

c 06.06.07 (BTD) changed interpretation of D2H
c                previously interpreted as ^2H/^1H
c                now interpreted as deuterium fraction ^2H/(^1H+^2H)
c                However: we are now only adding the C-D stretch
c                as other deuterated features are presumably already
c                included in the empirical representation of the
c                5-20um absorption cross section

c      FRACT_D=D2H/(1.D0+D2H)
      FRACT_D=D2H
c-----------------------------------------------------------------------

      FRACT_H=1.0-FRACT_D

      H2C=DBLE(NH)/DBLE(NC)

c set all enhancement factors here.
c 2006.09.09 (BTD) adjust all formulae to have enhancement factors = 1

      ENHANCE_033I=1.D0
      ENHANCE_033N=1.D0
      ENHANCE_052=1.D0
      ENHANCE_062=1.D0
c------------------------
c enhance 7.42, 7.60, 7.85um features by factor ENHANCE_077
      ENHANCE_077=1.D0

c------------------------
c enhance 8.33, 8.61um features by factor ENHANCE_086
c      ENHANCE_086=1.D0
c 080331 (BTD) v13:change ENHANCE_086 = 1 -> 1.1
c      ENHANCE_086=1.1D0
c 080408 (BTD) v14 change ENHANCE_086 = 1.1 -> 1.25
      ENHANCE_086=1.25D0
c------------------------
c enhance 10.68, 11.23, 11.30um features by factor ENHANCE_113
      ENHANCE_113=1.D0

      ENHANCE_120=1.0D0
c------------------------
c enhance 12.62, 12.69um features by factor ENHANCE_126
      ENHANCE_126=1.0D0
      ENHANCE_135=1.0D0
c 06.10.28 BTD:
c      ENHANCE_142=1.0D0
c      ENHANCE_159=1.0D0
c      ENHANCE_142=0.9D0
c      ENHANCE_159=0.8D0
c 06.10.31 BTD
      ENHANCE_142=1.0D0
      ENHANCE_159=1.0D0
      ENHANCE_164=1.0D0
c-------------------------
c enhance 17.0, 17.4, 17.9um features by factor ENHANCE_17
c 06.10.28 BTD:
c      ENHANCE_17=1.D0
c      ENHANCE_189=1.D0
c      ENHANCE_17=0.7D0
c      ENHANCE_189=0.7D0
c 06.10.31 BTD redefine
      ENHANCE_17=1.D0
      ENHANCE_189=1.D0
c 08.03.31 BTD v13 reduce 17um feature by 10%
c      ENHANCE_17=0.9D0
c 08.04.08 BTD v14 reduce 17um feature by another 5%
      ENHANCE_17=0.85D0
c------------------------
      WL_033=3.30D0
      WLD_033=3.30D0*CDSTRE_SHIFT

c FWHM: 0.20 um -- diffuse ISM (IRTS); 0.04 -- other sources;
c take gamma = .012 (fwhm=.012*3.3=.0396um)

      GAMMA_033=0.012D0

c 08.07.23 BTD v15: set FIP_033=1 for 3.3um C-H stretch.

      FIP_033=1.D0

      IF(STATE.EQ.'ionized')THEN
         INT_CROSSCT_033=89.4D-20*H2C*ENHANCE_033I*FRACT_H
         INT_CROSSCTD_033=89.4D-20*H2C*CD_REDUCED*FRACT_D*ENHANCE_033I
      ELSEIF(STATE.EQ.'neutral')THEN
         INT_CROSSCT_033=394.D-20*H2C*ENHANCE_033N*FRACT_H
         INT_CROSSCTD_033=394.D-20*H2C*CD_REDUCED*FRACT_D*ENHANCE_033N
      ENDIF

c add weak 5.27um feature as found in spectra of Orion Bar and M17
c (see SWS spectra in Peeters et al 2004, ApJ 604, 252,
c and SINGS galaxy spectra (Smith etal 2006)
c wavelength, gamma from Smith etal (2006)
c integrated cross section based on eyeball
c analysis of Orion Bar and M17 spectra, and SINGS spectr

c 2006.09.07: change wavelength from 5.25 to 5.27um
c             change gamma from 0.030 to 0.034

      WL_052=5.27D0
      GAMMA_052=0.034D0
      IF(STATE.EQ.'ionized')THEN
         INT_CROSSCT_052=20.D-20*ENHANCE_052
      ELSEIF(STATE.EQ.'neutral')THEN
         INT_CROSSCT_052=2.5D-20*ENHANCE_052
      ENDIF

c 2008.07.23 (BTD): 5.27um mode tentatively attributed to combination/overtones
c                   of out-of-plane vibrations (Allamandola et al 1989a)

      FIP_052=0.D0

c add weak 5.70um feature as found in SINGS spectra by Smith etal 2006:
c BTD eyeball estimate of integrated cross section scaled to 6.2um feature
c 2006.06.25 (BTD) v7: introduce ENHANCE_062 factor
c 2006.09.07 (BTD) v8: reduce gamma_057 from 0.040 to 0.035 (Smith etal 2006)
c                      reduce cross section based on sample of SINGS
c                      spectra which have median P(5.70)/P(6.22) approx 0.07
c                      for 12 galaxies with P(6.22)> 2e-7 W m-2 sr-1
c                      [but with large dispersion]

      WL_057=5.70D0
      GAMMA_057=0.035D0
      IF(STATE.EQ.'ionized')THEN
         INT_CROSSCT_057=32.D-20*ENHANCE_062
      ELSEIF(STATE.EQ.'neutral')THEN
         INT_CROSSCT_057=4.0D-20*ENHANCE_062
      ENDIF

c 2008.07.23 (BTD) 5.70um mode tentatively identified as combination/overtone
c                  modes of C_H out-of-plane bending modes
c                  (Allamandola et al 1989a)

      FIP_057=0.D0

      WL_062=6.22D0
      WLD_062=6.22D0*CC_SHIFT

c FWHM: 0.67 um for diffuse cloud (IRTS); 0.30 for other sources;
c take gamma = .032 (fwhm=.032*6.2=0.198)
c at some point changed this to 0.0284
c 2006.09.07 change gamma from .0284 to .030 [Smith etal 2006 Table 4]

      GAMMA_062=0.030D0

      IF(STATE.EQ.'ionized')THEN
         INT_CROSSCT_062=235.5D-20*ENHANCE_062*FRACT_H
         INT_CROSSCTD_062=235.5D-20*ENHANCE_062*FRACT_D
      ELSEIF(STATE .EQ. 'neutral')THEN
         INT_CROSSCT_062=29.4D-20*ENHANCE_062*FRACT_H
         INT_CROSSCTD_062=29.4D-20*ENHANCE_062*FRACT_D
      ENDIF

c 2008.07.28 (BTD) The 6.2um feature is assumed to be in-plane C-C stretch

      FIP_062=1.D0

c add weak 6.69um feature as found in SINGS spectra by Smith etal 2006:
c BTD eyeball estimate 0.2*integrated cross section scaled to 6.2um feature
c 2006.09.08 (BTD) 14 sings galaxies with p(6.22)> 2e-7 W m-2 sr-1 have
c                  median p(6.69)/p(6.22)=0.28
c                  guess this corresponds to cross section ratio 0.25

      WL_067=6.69D0
      GAMMA_067=0.070D0
      IF(STATE.EQ.'ionized')THEN
         INT_CROSSCT_067=0.25*INT_CROSSCT_062
      ELSEIF(STATE.EQ.'neutral')THEN
         INT_CROSSCT_067=0.25*INT_CROSSCT_062
      ENDIF

c 2008.07.28 (BTD) The 6.69um feature is unidentified (see Table 1 of DL07).  
c                  We will assume it to be an in-plane C-C vibration

      FIP_067=1.D0

c 05.10.31 based on analysis of SINGS spectra, we separate the 7.7um
c          feature into 3 components
c  wave   gamma  p/ptot
c 7.417   0.126  0.264 -> 0.26
c 7.598   0.044  0.390 -> 0.39
c 7.850   0.053  0.346 -> 0.35

c 06.09.08 analysis of 15 sings spectra with p(7.598)>2e-7 W m-2 sr-1
c  wave  p/ptot
c 7.417  1.15/(1.15+1+1.21)=0.342
c 7.598  1/(1.15+1+1.17)=   0.298
c 7.850  1.21/(15.1+1+1.21)=0.360

      WL_074=7.417D0
      WLD_074=WL_074*CC_SHIFT
      GAMMA_074=0.126D0
      F1=0.342D0

      WL_076=7.598D0
      WLD_076=7.598D0*CC_SHIFT
      GAMMA_076=0.044D0
      F2=0.298D0

      WL_079=7.850D0
      WLD_079=7.850D0*CC_SHIFT
      GAMMA_079=0.053
      F3=0.360D0

      IF(STATE .EQ.'ionized')THEN
         INT_CROSSCT_074=F1*548.D-20*ENHANCE_077*FRACT_H
         INT_CROSSCT_076=F2*548.D-20*ENHANCE_077*FRACT_H
         INT_CROSSCT_079=F3*548.D-20*ENHANCE_077*FRACT_H
         INT_CROSSCTD_074=F1*548.D-20*ENHANCE_077*FRACT_D
         INT_CROSSCTD_076=F2*548.D-20*ENHANCE_077*FRACT_D
         INT_CROSSCTD_079=F3*548.D-20*ENHANCE_077*FRACT_D
      ELSEIF(STATE .EQ. 'neutral')THEN
         INT_CROSSCT_074=F1*60.9D-20*ENHANCE_077*FRACT_H
         INT_CROSSCT_076=F2*60.9D-20*ENHANCE_077*FRACT_H
         INT_CROSSCT_079=F3*60.9D-20*ENHANCE_077*FRACT_H
         INT_CROSSCTD_074=F1*60.9D-20*ENHANCE_077*FRACT_D
         INT_CROSSCTD_076=F2*60.9D-20*ENHANCE_077*FRACT_D
         INT_CROSSCTD_079=F3*60.9D-20*ENHANCE_077*FRACT_D
      ENDIF

c 2008.07.23 (BTD) The 7.42, 7.60, 7.85 um features are assumed to be
c                  in-plane C-C stretch, mixed with in-plane C-H bending

      FIP_074=1.D0
      FIP_076=1.D0
      FIP_079=1.D0

c 05.11.28 based on analysis of SINGS spectrum of NGC 0337, 
c          we separate the 8.6um
c          feature into 2 components
c  wave   gamma p/ptot
c 8.33    .052  .208 -> .20
c 8.61    .039  .792 -> .80

      WL_083=8.33D0
      GAMMA_083=0.052D0
      F1=0.20D0
      WLD_083=WL_083*CDBEND_SHIFT

      WL_086=8.610D0
      GAMMA_086=0.039D0
      F2=0.80D0
      WLD_086=WL_086*CDBEND_SHIFT

c 2008.07.23 (BTD) 8.33 and 8.61um features are assumed to be mix of
c                  in-plane C-C stretch and in-plane C-H bending

      FIP_083=1.D0
      FIP_086=1.D0

      IF(STATE.EQ.'ionized')THEN
         INT_CROSSCT_083=F1*242.D-20*H2C*ENHANCE_086*FRACT_H
         INT_CROSSCTD_083=F1*242.D-20*H2C*ENHANCE_086*CD_REDUCED*FRACT_D
         INT_CROSSCT_086=F2*242.D-20*H2C*ENHANCE_086*FRACT_H
         INT_CROSSCTD_086=F2*242.D-20*H2C*ENHANCE_086*CD_REDUCED*FRACT_D
      ELSEIF(STATE.EQ.'neutral')THEN
         INT_CROSSCT_083=F1*34.7D-20*H2C*ENHANCE_086*FRACT_H
         INT_CROSSCTD_083=F1*34.7D-20*H2C*ENHANCE_086*CD_REDUCED*FRACT_D
         INT_CROSSCT_086=F2*34.7D-20*H2C*ENHANCE_086*FRACT_H
         INT_CROSSCTD_086=F2*34.7D-20*H2C*ENHANCE_086*CD_REDUCED*FRACT_D
      ENDIF

c 2006.09.07 (BTD) add new weak feature at 10.68um seen in SINGS spectra
c                  gamma=0.020
c                  for 10 galaxies with p(11.23)>3e-7 W m-2 sr-1
c                  median p(10.68)/p(11.23)=0.023
c                  guess sigma_int = 0.025*sigma_int(11.23)
c                                  = 0.3e-20*(H/C)

      WL_107=10.68D0
      GAMMA_107=0.020D0

c 2008.06.23 (BTD)
      IF(STATE.EQ.'ionized')THEN
         INT_CROSSCT_107=0.3D-20*H2C*ENHANCE_113*FRACT_H
      ELSEIF(STATE.EQ.'neutral')THEN
         INT_CROSSCT_107=0.3D-20*H2C*ENHANCE_113*FRACT_H
      ENDIF
         
c 2008.07.23 (BTD) 10.68um feature is assumed to be C-H out-of-plane bending
c                  (see Table 1 of DL07)

      FIP_107=0.D0

c we consider 3 C-H out-of-plane features: solo (11.3),
c duet (11.9), and trio (12.7). we assume equal C-H units for
c all three modes (suitable for large compact symmetric PAHs)
c (so there is a factor of 1/3).

c we separate the 11.3um feature into 2 components:
c wave  gamma  p/ptot
c 11.23 .010   0.182 -> 0.18
c 11.29 .029   0.818 -> 0.82

c 2006.09.07 (BTD) change gamma_112 from 0.010 to 0.012 (Smith etal 2006)
c 2006.09.08 (BTD) 15 sings galaxies with p(11.33>2e-7 W m-2 sr-1) have
c                  median p(11.23)/p(11.33)=0.36
c                  therefore change f1 to 0.36/(1+0.36)=0.265
c                                   f2 to 1/(1+0.36)=0.735

      WL_112=11.23D0
      WLD_112=11.23D0*CDBEND_SHIFT
      GAMMA_112=0.012
      F1=0.265D0

c 2006.09.07 (BTD) change wl_113 from 11.30 to 11.33
c                  change gamma_113 from 0.029 to 0.032
c 2006.11.01 (BTD) change wl_113 from 11.30 to 11.33 (didn't do it before)
      WL_113=11.33D0
      WLD_113=WL_113*CDBEND_SHIFT
      GAMMA_113=0.032
      F2=0.735D0

      IF(STATE.EQ.'ionized')THEN
         INT_CROSSCT_112=F1*66.67D-20*ENHANCE_113*H2C*FRACT_H
         INT_CROSSCTD_112=F1*66.67D-20*ENHANCE_113*H2C*
     &                    CD_REDUCED*FRACT_D
         INT_CROSSCT_113=F2*66.67D-20*ENHANCE_113*H2C*FRACT_H
         INT_CROSSCTD_113=F2*66.67D-20*ENHANCE_113*H2C*
     &                    CD_REDUCED*FRACT_D
      ELSEIF(STATE.EQ.'neutral')THEN
         INT_CROSSCT_112=F1*71.2D-20*ENHANCE_113*H2C*FRACT_H
         INT_CROSSCTD_112=F1*71.2D-20*ENHANCE_113*H2C*
     &                    CD_REDUCED*FRACT_D
         INT_CROSSCT_113=F2*71.2D-20*ENHANCE_113*H2C*FRACT_H
         INT_CROSSCTD_113=F2*71.2D-20*ENHANCE_113*H2C*
     &                    CD_REDUCED*FRACT_D
      ENDIF

c 2008.07.23 (BTD) 11.23, 11.33 features are assumed to
c                  be C-H solo out-of-plane bending (see DL07 Table 1)

      FIP_112=0.D0
      FIP_113=0.D0

c wavelength=11.99um and gamma=0.050 from SINGS data
c integrated cross sections from Li & Draine 2001, Table 1
c 2006.09.07 (BTD) change gamma from 0.050 to 0.045 (Smith etal 2006)

      WL_120=11.99D0
      WLD_120=WL_120*CDBEND_SHIFT
      GAMMA_120=0.045

      IF(STATE.EQ.'ionized')THEN
         INT_CROSSCT_120=20.5D-20*ENHANCE_120*H2C*FRACT_H
         INT_CROSSCTD_120=20.5D-20*ENHANCE_120*H2C*CD_REDUCED*FRACT_D
      ELSEIF(STATE.EQ.'neutral')THEN
         INT_CROSSCT_120=24.2D-20*ENHANCE_120*H2C*FRACT_H
         INT_CROSSCTD_120=24.2D-20*ENHANCE_120*H2C*CD_REDUCED*FRACT_D
      ENDIF

c 2008.07.23 (BTD) 11.99um feature is asssumed to be duo C-H out-of-plane
c                  bending mode (see DL07 Table 1)

      FIP_120=0.D0

c wavelength=12.61um and gamma=0.0435 from SINGS data
c integrated cross sections from Li & Draine 2001, Table 1
c multiplied by enhancement factor ENHANCE_126 (set above to 1.25)
c to increase strength of 12.6um feature to agree with SINGS spectra
c 2006.09.07: reset ENHANCE_126 to 1.0 with corresponding change in
c             equation to calculate INT_CROSSCT_126
c             change wavelength from 12.61 to 12.62
c             change gamma from .0435 to .042

      WL_126=12.62D0
      WLD_126=WL_126*CDBEND_SHIFT
      GAMMA_126=.042D0

      IF(STATE.EQ.'ionized')THEN
         INT_CROSSCT_126=ENHANCE_126*31.0D-20*H2C*FRACT_H
         INT_CROSSCTD_126=ENHANCE_126*31.0D-20*H2C*CD_REDUCED*FRACT_D
      ELSEIF(STATE.EQ.'neutral')THEN
         INT_CROSSCT_126=ENHANCE_126*34.8D-20*H2C*FRACT_H
         INT_CROSSCTD_126=ENHANCE_126*34.8D-20*H2C*CD_REDUCED*FRACT_D
      ENDIF

c 2008.07.23 (BTD) 12.62um feature is assumed to be trio C-H out-of-plane
c                  bending mode (see DL07 Table 1)

      FIP_126=0.D0

c 2006.09.07: add new narrow feature at 12.69um from Smith etal 2006
c             gamma=0.013
c             16 sings galaxies with P(12.62) > 1e-7
c             median p(12.69)/p(12.62)=0.04
c             estimate sigma_int(12.69)=0.04*sigma_int(12.62)
c             with sigma_int(12.62)=34.79d-20(H/C) for neutrals, 
c                                   31.04d-20(H/C) for neutrals
c             we simply adopt 1.3d-20(H/C) for both
      WL_127=12.69D0
      WLD_127=WL_127*CDBEND_SHIFT
      GAMMA_127=0.013D0
      INT_CROSSCT_127=ENHANCE_126*1.3D-20*H2C*FRACT_H
      INT_CROSSCTD_127=ENHANCE_126*1.3D-20*H2C*CD_REDUCED*FRACT_D

c 2008.07.23 (BTD) 12.69um feature is considered to be trio out-of-plane
c                  bending mode (DL07 Table 1)

      FIP_127=0.D0

c wavelength=13.60um and gamma=0.020 from SINGS data
c integrated cross section is a guess...
c 2006.09.07: change wavelength from 13.60 to 13.48
c             change gamma from 0.020 to 0.040
c             relate cross section to cross section of 12.62um feature
c             based on 0.30 = power in 13.48um feature/
c                             power in 12.62um feature
c             adopt same cross section for neutral and ion

      WL_135=13.48D0
      GAMMA_135=0.040D0
      INT_CROSSCT_135=8.D-20*ENHANCE_135*H2C

c 2008.07.23 (BTD) 13.48um feature is tentativeley assumed to be quartet
c                  out-of-plane C-H bending mode (DL07 Table 1)

      FIP_135=0.D0

c for wavelengths > 14um, we do not consider effects of D, since
c features are presumably C-C bending modes

c wavelength=14.19um and gamma=0.025 from SINGS data
c 06.06.07 (BTD) changed from 0.5d-20 -> 1.0d-20 
c 06.09.07 power(14.19um)/power(12.62um) = 0.12
c          for H/C=0.25, sigmaint(12.62um)=15.5e-20 cm/C for ion
c          therefore estimate sigmaint= 1.86e-20 cm/C for ion
c                              adopt same value for neutral
c 08.03.17 (BTD) v12: increase strength of 14.19um feature 
c          from 0.45e-20 to 0.60e-20

      WL_142=14.19D0
      GAMMA_142=0.025D0
      INT_CROSSCT_142=0.60D-20*ENHANCE_142

c 2008.07.23 (BTD) 14.19um feature is tentatively assigned to quartet
c                  out-of-plane bending mode (see DL07 Table 1)

      FIP_142=0.D0

c weak feature at 15.90um seen in SINGS spectra (Smith etal 2006)
c BTD eyeball estimate for integrated cross section
c 2006.09.07 (BTD) changed estimate for sigmaint(15.9)
c            previous estimate was sigmaint(15.9)=0.75d-20 cm/C
c            but for 13 sings galaxies with P(16.45)>1e-8 W m-2 sr-1
c            median P(15.90)/P(16.45)=0.10 (with considerable scatter)
c            sigma_int(16.45)=0.75e-20 cm/C
c            therefore sigma_int(15.90)=0.10*0.75e-20=0.08e-20 cm/C

      WL_159=15.90D0
      GAMMA_159=0.020D0
      INT_CROSSCT_159=0.04D-20*ENHANCE_159

c 2008.07.23 (BTD) 15.90um feature is unidentified.  
c                  We assume 2/3 in-plane vibrations, 1/3 to out-of-plane

      FIP_159=(2.D0/3.D0)

c wavelength=16.447um and gamma=0.014 from SINGS data
c integrated cross section 5.52d-20cm taken from 16.4um feature in
c Moutou et al. 1996 and Moutou et al. 2000: 16.4 micron, FWHM=0.16um.

      WL_164=16.447
      GAMMA_164=0.014

c we reduce the integrated cross section from 5.52D-20 in LD2001
c                                        to   3.00D-20
c      INT_CROSSCT_164=5.52D-20
c      INT_CROSSCT_164=3.00D-20
c 05.12.05 (BTD) further reduce to 2.00D-20
c 06.06.25 (BTD) further reduce to 0.75d-20 (enhance_17=0.5)
c 06.09.06 (BTD) keep at 0.75d-20 for enhance_17=1
c 06.10.31 (BTD) reset to 0.5d-20 for enhance_164=1

      INT_CROSSCT_164=0.5D-20*ENHANCE_164

c 2008.07.23 (BTD) 16.45um feature is tentatively assigned to C-C-C bending 
c                  modes.  The character of these vibrations is uncertain
c                  Tentatively take them to be 2/3 in-plane, 1/3 out-of-plane,
c                  but need to see whether future theoretical calculations 
c                  support this.

      FIP_164=(2.D0/3.D0)

c three components to 17um complex
c  wave  gamma p/ptot
c 17.038 0.065  0.504
c 17.377 0.012  0.247
c 17.873 0.016  0.250

c      F1=.50D0
c      F2=.25D0
c      F3=.25D0

c above was in error.  From NGC0337 and NGC6946 we take
c      F1=0.88
c      F2=0.08
c      F3=0.04
c
c 2006.09.07 from analysis of PAHFit results for 31 galaxies:
c      F1=1/(1+.05+.03)=  0.926
c      F2=.05/(1+.05+.03)=0.046
c      F3=.03/(1+.05+.03)=0.028

      F1=0.926
      F2=0.046
      F3=0.028

c following integrated cross section (summed over 3 components) is a guess.
c 
c      INT_CROSSCT_17=2.4D-20*ENHANCE_17
c 080104 (BTD) v10: try larger guess
c      INT_CROSSCT_17=3.5D-20*ENHANCE_17
c 080117 (BTD) v11: try larger guess
c      INT_CROSSCT_17=4.5D-20*ENHANCE_17
c 080317 (BTD) v12: try smaller guess

      INT_CROSSCT_17=3.8D-20*ENHANCE_17
      WL_170=17.04D0
      GAMMA_170=0.065D0
      INT_CROSSCT_170=F1*INT_CROSSCT_17

      WL_174=17.375D0
      GAMMA_174=0.012D0
      INT_CROSSCT_174=F2*INT_CROSSCT_17

      WL_179=17.87D0
      GAMMA_179=0.016D0
      INT_CROSSCT_179=F3*INT_CROSSCT_17

c 2008.07.23 (BTD) 17.04, 17.375, and 17.87um features are tentatively
c                  asssigned to C-C-C bending modes.  The character of these
c                  modes is uncertain.
c                  Tentatively take them to be 2/3 in-plane, 1/3 out-of-plane,
c                  but need to see whether future theoretical calculations 
c                  support this.

      FIP_170=(2.D0/3.D0)
      FIP_174=(2.D0/3.D0)
      FIP_179=(2.D0/3.D0)

c 2006.09.07 add new feature at 18.92um
c            15 SINGS galaxies with P(17.04) > 1e-7 W cm-2 sr-1
c            have median P(18.92um)/P(17.04)=0.047
c            P(18.92)/P(17um complex)=.047*0.926=.044
c            Therefore try 
c            sigma_int(18.92)=0.04*sigma_int(17um complex)
c                            =0.14d-20*ENHANCE_189

      WL_189=18.92D0
      GAMMA_189=0.019D0
      INT_CROSSCT_189=0.10D-20*ENHANCE_189

c 2008.07.23 (BTD) 18.92um feature is tentatively assigned to C-C-C bending
c                  modes.  The character of these vibrations is uncertain
c                  Tentatively take them to be 2/3 in-plane, 1/3 out-of-plane,
c                  but need to see whether future theoretical calculations 
c                  support this

      FIP_189=(2.D0/3.D0)

c Add a FIR continuum (the tail of Drude profile).
c Drude profile (for FIR): gamma=FWHM, wlpeak;
c Li & Draine 2001 took wl=26um, gamma=0.7
c Here we experiment with wl=20um, gamma=0.7

c experiment to reproduce SINGS photometry:

c      WL_FIR=17.0D0
c      GAMMA_FIR=0.7D0
c      INT_CROSSCT_FIR=50.D-20

c 051205 (BTD) experiment with further increase in FIR continuum
c              seeking to fill in trough between 12.6um and 16.4um
c 060415 (BTD) new experiment to raise floor of valley between
c              12.7 and 16.4 : 
c              shift WL_FIR to 15um
c              broaden by changing gamma from 0.7 to 0.8
c 060415 (BTD) add ENHANCE_17 enhancement factor

      WL_FIR=15.0D0
      GAMMA_FIR=0.8D0

c      INT_CROSSCT_FIR=70.D-20*ENHANCE_17

c 060418 (BTD) raise INT_CROSSCT_FIR from 1.4e-18 to 2.0d-18

c      INT_CROSSCT_FIR=50.D-20*ENHANCE_17
c 080104 (BTD) v10 lower INT_CROSSSCT_FIR from 50.d-20 to 30.d-20
c      INT_CROSSCT_FIR=30.D-20*ENHANCE_17
c 080117 (BTD) v11 lower INT_CROSSSCT_FIR from 30.d-20 to 25.d-20
      INT_CROSSCT_FIR=25.D-20*ENHANCE_17

c 080723 (BTD) v15 FIR emission is likely to involve both in-plane
c                  and out-of-plane modes.
c                  Tentatively take these to be 2/3 in-plane,
c                  1/3 out-of-plane

      FIP_FIR=(2.D0/3.D0)

      DCRS=INT_CROSSCT_033*DRUDES(WL_033,GAMMA_033,WAVEUM)
      CRS_SCT_FEATURE=DCRS
      CRS_IP=CRS_IP+FIP_033*DCRS

      DCRS=INT_CROSSCT_052*DRUDES(WL_052,GAMMA_052,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_052*DCRS

      DCRS=INT_CROSSCT_057*DRUDES(WL_057,GAMMA_057,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_057*DCRS

      DCRS=INT_CROSSCT_062*DRUDES(WL_062,GAMMA_062,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_062*DCRS

      DCRS=INT_CROSSCT_067*DRUDES(WL_067,GAMMA_067,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_067*DCRS

      DCRS=INT_CROSSCT_074*DRUDES(WL_074,GAMMA_074,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_074*DCRS

      DCRS=INT_CROSSCT_076*DRUDES(WL_076,GAMMA_076,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_076*DCRS

      DCRS=INT_CROSSCT_079*DRUDES(WL_079,GAMMA_079,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_079*DCRS

      DCRS=INT_CROSSCT_083*DRUDES(WL_083,GAMMA_083,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_083*DCRS

      DCRS=INT_CROSSCT_086*DRUDES(WL_086,GAMMA_086,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_086*DCRS

      DCRS=INT_CROSSCT_107*DRUDES(WL_107,GAMMA_107,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_107*DCRS

      DCRS=INT_CROSSCT_112*DRUDES(WL_112,GAMMA_112,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_112*DCRS

      DCRS=INT_CROSSCT_113*DRUDES(WL_113,GAMMA_113,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_113*DCRS

      DCRS=INT_CROSSCT_120*DRUDES(WL_120,GAMMA_120,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_120*DCRS

      DCRS=INT_CROSSCT_126*DRUDES(WL_126,GAMMA_126,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_126*DCRS

      DCRS=INT_CROSSCT_127*DRUDES(WL_127,GAMMA_127,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_127*DCRS

      DCRS=INT_CROSSCT_135*DRUDES(WL_135,GAMMA_135,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_135*DCRS

      DCRS=INT_CROSSCT_142*DRUDES(WL_142,GAMMA_142,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_142*DCRS

      DCRS=INT_CROSSCT_159*DRUDES(WL_159,GAMMA_159,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_159*DCRS

      DCRS=INT_CROSSCT_164*DRUDES(WL_164,GAMMA_164,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_164*DCRS

      DCRS=INT_CROSSCT_170*DRUDES(WL_170,GAMMA_170,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_170*DCRS

      DCRS=INT_CROSSCT_174*DRUDES(WL_174,GAMMA_174,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_174*DCRS

      DCRS=INT_CROSSCT_179*DRUDES(WL_179,GAMMA_179,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_179*DCRS

      DCRS=INT_CROSSCT_189*DRUDES(WL_189,GAMMA_189,WAVEUM)
      CRS_SCT_FEATURE=CRS_SCT_FEATURE+DCRS
      CRS_IP=CRS_IP+FIP_189*DCRS

      DCRS=INT_CROSSCT_FIR*DRUDES(WL_FIR,GAMMA_FIR,WAVEUM)
      CRS_SCT_FIRCONT=DCRS
      CRS_IP=CRS_IP+FIP_FIR*DCRS

      IF(D2H.GT.0.D0)THEN
         CRS_SCT_FEATURE=CRS_SCT_FEATURE+
     &             INT_CROSSCTD_033*DRUDES(WLD_033,GAMMA_033,WAVEUM)
c     &             +INT_CROSSCTD_062*DRUDES(WLD_062,GAMMA_062,WAVEUM)+
c     &             INT_CROSSCTD_074*DRUDES(WLD_074,GAMMA_074,WAVEUM)+
c     &             INT_CROSSCTD_076*DRUDES(WLD_076,GAMMA_076,WAVEUM)+
c     &             INT_CROSSCTD_079*DRUDES(WLD_079,GAMMA_079,WAVEUM)+
c     &             INT_CROSSCTD_083*DRUDES(WLD_083,GAMMA_083,WAVEUM)+
c     &             INT_CROSSCTD_086*DRUDES(WLD_086,GAMMA_086,WAVEUM)+
c     &             INT_CROSSCTD_112*DRUDES(WLD_112,GAMMA_112,WAVEUM)+
c     &             INT_CROSSCTD_113*DRUDES(WLD_113,GAMMA_113,WAVEUM)+
c     &             INT_CROSSCTD_120*DRUDES(WLD_120,GAMMA_120,WAVEUM)+
c     &             INT_CROSSCTD_126*DRUDES(WLD_126,GAMMA_126,WAVEUM)+
C     &             INT_CROSSCTD_127*DRUDES(WLD_127,GAMMA_127,WAVEUM)
      ENDIF
      CROSS_SECTION=CRS_SCT_UV+CRS_SCT_NIR+CRS_SCT_FEATURE+
     &              CRS_SCT_FIRCONT
      FIP=CRS_IP/CROSS_SECTION

c cross_section: cm^2/C-atom.
c*** diagnostic
c      if(wavelengthcm.gt.2.1e-5.and.wavelengthcm.lt.2.2e-5)then
c         write(0,*)'return from pah_crs_sct_v9, CROSS_SECTION=',
c     &             CROSS_SECTION
c      endif
c***
      RETURN
      END

      FUNCTION DRUDES(WAVE0,GAMMA,WAVE)
      IMPLICIT NONE
      DOUBLE PRECISION DRUDES,GAMMA,WAVE,WAVE0

c given
c     WAVE0 = resonance wavelength (um)
c     GAMMA = dimensionless broadening parameter
c     WAVE  = wavelength (um)
c
c returns
c     DRUDES = (2/pi)*gamma*wave0/((wave/wave0-wave0/wave)^2+gamma^2
c               (cm)
c
c 0.6366197723D-4 = (2/pi)*micron/cm
c
c B.T. Draine, Osservatorio Astrofisico di Arcetri
c 2005.10.31
c-----------------------------------------------------------------------

      DRUDES=0.6366197723D-4*GAMMA*WAVE0/
     &       (GAMMA**2+(WAVE/WAVE0-WAVE0/WAVE)**2)
      RETURN
      END
      SUBROUTINE REFICE(IUNIT,XLAM,T,RN,CN,ABSIND,ABSCOF)
C Arguments:
      INTEGER IUNIT
      REAL ABSCOF,ABSIND,CN,XLAM,RN,T
C Parameters:
      INTEGER I,LT1,LT2,NWL,NWLT
      PARAMETER(NWL=468,NWLT=62)
C Local variables:
      REAL
     &   ALAM,CUTICE,PI,T1,T2,TK,WLMAX,WLMIN,
     &   X,X1,X2,Y,Y1,Y2,YLO,YHI

      REAL
     &   TABIM(NWL),TABIMT(NWLT,4),TABRE(NWL),TABRET(NWLT,4),TEMREF(4),
     &   WL(NWL),WLT(NWLT)
C
C     DEFINES WAVELENGTH DEPENDENT COMPLEX INDEX OF REFRACTION FOR ICE.
C     ALLOWABLE WAVELENGTH RANGE EXTENDS FROM 0.045 MICRONS TO 8.6 METER
C     TEMPERATURE DEPENDENCE ONLY CONSIDERED BEYOND 167 MICRONS.
C
C     INTERPOLATION IS DONE     RN  VS. LOG(XLAM)
C                               RN  VS.        T
C                           LOG(CN) VS. LOG(XLAM)
C                           LOG(CN) VS.        T
C
C     STEPHEN G. WARREN - 1983
C     DEPT. OF ATMOSPHERIC SCIENCES
C     UNIVERSITY OF WASHINGTON
C     SEATTLE, WA  98195
C
C     BASED ON
C
C        WARREN,S.G.,1984.
C        OPTICAL CONSTANTS OF ICE FROM THE ULTRAVIOLET TO THE MICROWAVE.
C        APPLIED OPTICS,23,1206-1225
C
C     INPUT PARAMETERS
C
C     IUNIT = 0 FOR WAVELENGTH SPECIFIED IN MICRONS
C           = 1 FOR WAVELENGTH SPECIFIED IN MILLIMETERS
C           = 2 FOR WAVELENGTH SPECIFIED IN CENTIMETERS
C           = 3 FOR WAVELENGTH SPECIFIED IN INVERSE CENTIMETERS ( WAVE N
C     XLAM = WAVELENGTH ( MICRONS OR MM OR CM OR CM**-1 )
C     T = TEMPERATURE ( DEGREES KELVIN )
C
C     OUTPUT PARAMETERS
C
C     RN = REAL PORTION ( SCATTERING )
C     CN = COMPLEX PORTION ( ABSORPTION )
C     ABSIND = ABSORPTIVE INDEX ( CN/RN )
C     ABSCOF = ABORPTION COEFFICIENT ( 4*PI*CN/XLAM )
C
C      DIMENSION WL(NWL),WLT(NWLT)
C      DIMENSION TABRE(NWL),TABRET(NWLT,4),TABIM(NWL),TABIMT(NWLT,4)
C      DIMENSION TEMREF(4)
C
C     REFERENCE TEMPERATURES ARE -1.0,-5.0,-20.0, AND -60.0 DEG CENTIGRA
C
      DATA TEMREF/272.16,268.16,253.16,213.16/
C
      DATA WLMIN,WLMAX/0.045,8.6E6/
      DATA CUTICE/167.0/
C
      DATA (WL(I),I=1,114)/
     +0.4430E-01,0.4510E-01,0.4590E-01,0.4680E-01,0.4770E-01,0.4860E-01,
     +0.4960E-01,0.5060E-01,0.5170E-01,0.5280E-01,0.5390E-01,0.5510E-01,
     +0.5640E-01,0.5770E-01,0.5900E-01,0.6050E-01,0.6200E-01,0.6360E-01,
     +0.6530E-01,0.6700E-01,0.6890E-01,0.7080E-01,0.7290E-01,0.7380E-01,
     +0.7510E-01,0.7750E-01,0.8000E-01,0.8270E-01,0.8550E-01,0.8860E-01,
     +0.9180E-01,0.9300E-01,0.9540E-01,0.9920E-01,0.1033E+00,0.1078E+00,
     +0.1100E+00,0.1127E+00,0.1140E+00,0.1181E+00,0.1210E+00,0.1240E+00,
     +0.1272E+00,0.1295E+00,0.1305E+00,0.1319E+00,0.1333E+00,0.1348E+00,
     +0.1362E+00,0.1370E+00,0.1378E+00,0.1387E+00,0.1393E+00,0.1409E+00,
     +0.1425E+00,0.1435E+00,0.1442E+00,0.1450E+00,0.1459E+00,0.1468E+00,
     +0.1476E+00,0.1480E+00,0.1485E+00,0.1494E+00,0.1512E+00,0.1531E+00,
     +0.1540E+00,0.1550E+00,0.1569E+00,0.1580E+00,0.1589E+00,0.1610E+00,
     +0.1625E+00,0.1648E+00,0.1669E+00,0.1692E+00,0.1713E+00,0.1737E+00,
     +0.1757E+00,0.1779E+00,0.1802E+00,0.1809E+00,0.1821E+00,0.1833E+00,
     +0.1843E+00,0.1850E+00,0.1860E+00,0.1870E+00,0.1880E+00,0.1890E+00,
     +0.1900E+00,0.1910E+00,0.1930E+00,0.1950E+00,0.2100E+00,0.2500E+00,
     +0.3000E+00,0.3500E+00,0.4000E+00,0.4100E+00,0.4200E+00,0.4300E+00,
     +0.4400E+00,0.4500E+00,0.4600E+00,0.4700E+00,0.4800E+00,0.4900E+00,
     +0.5000E+00,0.5100E+00,0.5200E+00,0.5300E+00,0.5400E+00,0.5500E+00/
      DATA (WL(I),I=115,228)/
     +0.5600E+00,0.5700E+00,0.5800E+00,0.5900E+00,0.6000E+00,0.6100E+00,
     +0.6200E+00,0.6300E+00,0.6400E+00,0.6500E+00,0.6600E+00,0.6700E+00,
     +0.6800E+00,0.6900E+00,0.7000E+00,0.7100E+00,0.7200E+00,0.7300E+00,
     +0.7400E+00,0.7500E+00,0.7600E+00,0.7700E+00,0.7800E+00,0.7900E+00,
     +0.8000E+00,0.8100E+00,0.8200E+00,0.8300E+00,0.8400E+00,0.8500E+00,
     +0.8600E+00,0.8700E+00,0.8800E+00,0.8900E+00,0.9000E+00,0.9100E+00,
     +0.9200E+00,0.9300E+00,0.9400E+00,0.9500E+00,0.9600E+00,0.9700E+00,
     +0.9800E+00,0.9900E+00,0.1000E+01,0.1010E+01,0.1020E+01,0.1030E+01,
     +0.1040E+01,0.1050E+01,0.1060E+01,0.1070E+01,0.1080E+01,0.1090E+01,
     +0.1100E+01,0.1110E+01,0.1120E+01,0.1130E+01,0.1140E+01,0.1150E+01,
     +0.1160E+01,0.1170E+01,0.1180E+01,0.1190E+01,0.1200E+01,0.1210E+01,
     +0.1220E+01,0.1230E+01,0.1240E+01,0.1250E+01,0.1260E+01,0.1270E+01,
     +0.1280E+01,0.1290E+01,0.1300E+01,0.1310E+01,0.1320E+01,0.1330E+01,
     +0.1340E+01,0.1350E+01,0.1360E+01,0.1370E+01,0.1380E+01,0.1390E+01,
     +0.1400E+01,0.1410E+01,0.1420E+01,0.1430E+01,0.1440E+01,0.1449E+01,
     +0.1460E+01,0.1471E+01,0.1481E+01,0.1493E+01,0.1504E+01,0.1515E+01,
     +0.1527E+01,0.1538E+01,0.1563E+01,0.1587E+01,0.1613E+01,0.1650E+01,
     +0.1680E+01,0.1700E+01,0.1730E+01,0.1760E+01,0.1800E+01,0.1830E+01,
     +0.1840E+01,0.1850E+01,0.1855E+01,0.1860E+01,0.1870E+01,0.1890E+01/
      DATA (WL(I),I=229,342)/
     +0.1905E+01,0.1923E+01,0.1942E+01,0.1961E+01,0.1980E+01,0.2000E+01,
     +0.2020E+01,0.2041E+01,0.2062E+01,0.2083E+01,0.2105E+01,0.2130E+01,
     +0.2150E+01,0.2170E+01,0.2190E+01,0.2220E+01,0.2240E+01,0.2245E+01,
     +0.2250E+01,0.2260E+01,0.2270E+01,0.2290E+01,0.2310E+01,0.2330E+01,
     +0.2350E+01,0.2370E+01,0.2390E+01,0.2410E+01,0.2430E+01,0.2460E+01,
     +0.2500E+01,0.2520E+01,0.2550E+01,0.2565E+01,0.2580E+01,0.2590E+01,
     +0.2600E+01,0.2620E+01,0.2675E+01,0.2725E+01,0.2778E+01,0.2817E+01,
     +0.2833E+01,0.2849E+01,0.2865E+01,0.2882E+01,0.2899E+01,0.2915E+01,
     +0.2933E+01,0.2950E+01,0.2967E+01,0.2985E+01,0.3003E+01,0.3021E+01,
     +0.3040E+01,0.3058E+01,0.3077E+01,0.3096E+01,0.3115E+01,0.3135E+01,
     +0.3155E+01,0.3175E+01,0.3195E+01,0.3215E+01,0.3236E+01,0.3257E+01,
     +0.3279E+01,0.3300E+01,0.3322E+01,0.3345E+01,0.3367E+01,0.3390E+01,
     +0.3413E+01,0.3436E+01,0.3460E+01,0.3484E+01,0.3509E+01,0.3534E+01,
     +0.3559E+01,0.3624E+01,0.3732E+01,0.3775E+01,0.3847E+01,0.3969E+01,
     +0.4099E+01,0.4239E+01,0.4348E+01,0.4387E+01,0.4444E+01,0.4505E+01,
     +0.4547E+01,0.4560E+01,0.4580E+01,0.4719E+01,0.4904E+01,0.5000E+01,
     +0.5100E+01,0.5200E+01,0.5263E+01,0.5400E+01,0.5556E+01,0.5714E+01,
     +0.5747E+01,0.5780E+01,0.5814E+01,0.5848E+01,0.5882E+01,0.6061E+01,
     +0.6135E+01,0.6250E+01,0.6289E+01,0.6329E+01,0.6369E+01,0.6410E+01/
      DATA (WL(I),I=343,456)/
     +0.6452E+01,0.6494E+01,0.6579E+01,0.6667E+01,0.6757E+01,0.6897E+01,
     +0.7042E+01,0.7143E+01,0.7246E+01,0.7353E+01,0.7463E+01,0.7576E+01,
     +0.7692E+01,0.7812E+01,0.7937E+01,0.8065E+01,0.8197E+01,0.8333E+01,
     +0.8475E+01,0.8696E+01,0.8929E+01,0.9091E+01,0.9259E+01,0.9524E+01,
     +0.9804E+01,0.1000E+02,0.1020E+02,0.1031E+02,0.1042E+02,0.1053E+02,
     +0.1064E+02,0.1075E+02,0.1087E+02,0.1100E+02,0.1111E+02,0.1136E+02,
     +0.1163E+02,0.1190E+02,0.1220E+02,0.1250E+02,0.1282E+02,0.1299E+02,
     +0.1316E+02,0.1333E+02,0.1351E+02,0.1370E+02,0.1389E+02,0.1408E+02,
     +0.1429E+02,0.1471E+02,0.1515E+02,0.1538E+02,0.1563E+02,0.1613E+02,
     +0.1639E+02,0.1667E+02,0.1695E+02,0.1724E+02,0.1818E+02,0.1887E+02,
     +0.1923E+02,0.1961E+02,0.2000E+02,0.2041E+02,0.2083E+02,0.2222E+02,
     +0.2260E+02,0.2305E+02,0.2360E+02,0.2460E+02,0.2500E+02,0.2600E+02,
     +0.2857E+02,0.3100E+02,0.3333E+02,0.3448E+02,0.3564E+02,0.3700E+02,
     +0.3824E+02,0.3960E+02,0.4114E+02,0.4276E+02,0.4358E+02,0.4458E+02,
     +0.4550E+02,0.4615E+02,0.4671E+02,0.4736E+02,0.4800E+02,0.4878E+02,
     +0.5003E+02,0.5128E+02,0.5275E+02,0.5350E+02,0.5424E+02,0.5500E+02,
     +0.5574E+02,0.5640E+02,0.5700E+02,0.5746E+02,0.5840E+02,0.5929E+02,
     +0.6000E+02,0.6100E+02,0.6125E+02,0.6250E+02,0.6378E+02,0.6467E+02,
     +0.6558E+02,0.6655E+02,0.6760E+02,0.6900E+02,0.7053E+02,0.7300E+02/
      DATA (WL(I),I=457,468)/
     +0.7500E+02,0.7629E+02,0.8000E+02,0.8297E+02,0.8500E+02,0.8680E+02,
     +0.9080E+02,0.9517E+02,0.1000E+03,0.1200E+03,0.1500E+03,0.1670E+03/
      DATA  WLT/
     +                                 0.1670E+03,0.1778E+03,0.1884E+03,
     +0.1995E+03,0.2113E+03,0.2239E+03,0.2371E+03,0.2512E+03,0.2661E+03,
     +0.2818E+03,0.2985E+03,0.3162E+03,0.3548E+03,0.3981E+03,0.4467E+03,
     +0.5012E+03,0.5623E+03,0.6310E+03,0.7943E+03,0.1000E+04,0.1259E+04,
     +0.2500E+04,0.5000E+04,0.1000E+05,0.2000E+05,0.3200E+05,0.3500E+05,
     +0.4000E+05,0.4500E+05,0.5000E+05,0.6000E+05,0.7000E+05,0.9000E+05,
     +0.1110E+06,0.1200E+06,0.1300E+06,0.1400E+06,0.1500E+06,0.1600E+06,
     +0.1700E+06,0.1800E+06,0.2000E+06,0.2500E+06,0.2900E+06,0.3200E+06,
     +0.3500E+06,0.3800E+06,0.4000E+06,0.4500E+06,0.5000E+06,0.6000E+06,
     +0.6400E+06,0.6800E+06,0.7200E+06,0.7600E+06,0.8000E+06,0.8400E+06,
     +0.9000E+06,0.1000E+07,0.2000E+07,0.5000E+07,0.8600E+07/
      DATA (TABRE(I),I=1,114)/
     +   0.83441,   0.83676,   0.83729,   0.83771,   0.83827,   0.84038,
     +   0.84719,   0.85522,   0.86047,   0.86248,   0.86157,   0.86093,
     +   0.86419,   0.86916,   0.87764,   0.89296,   0.91041,   0.93089,
     +   0.95373,   0.98188,   1.02334,   1.06735,   1.11197,   1.13134,
     +   1.15747,   1.20045,   1.23840,   1.27325,   1.32157,   1.38958,
     +   1.41644,   1.40906,   1.40063,   1.40169,   1.40934,   1.40221,
     +   1.39240,   1.38424,   1.38075,   1.38186,   1.39634,   1.40918,
     +   1.40256,   1.38013,   1.36303,   1.34144,   1.32377,   1.30605,
     +   1.29054,   1.28890,   1.28931,   1.30190,   1.32025,   1.36302,
     +   1.41872,   1.45834,   1.49028,   1.52128,   1.55376,   1.57782,
     +   1.59636,   1.60652,   1.61172,   1.61919,   1.62522,   1.63404,
     +   1.63689,   1.63833,   1.63720,   1.63233,   1.62222,   1.58269,
     +   1.55635,   1.52453,   1.50320,   1.48498,   1.47226,   1.45991,
     +   1.45115,   1.44272,   1.43498,   1.43280,   1.42924,   1.42602,
     +   1.42323,   1.42143,   1.41897,   1.41660,   1.41434,   1.41216,
     +   1.41006,   1.40805,   1.40423,   1.40067,   1.38004,   1.35085,
     +   1.33394,   1.32492,   1.31940,   1.31854,   1.31775,   1.31702,
     +   1.31633,   1.31569,   1.31509,   1.31452,   1.31399,   1.31349,
     +   1.31302,   1.31257,   1.31215,   1.31175,   1.31136,   1.31099/
      DATA (TABRE(I),I=115,228)/
     +   1.31064,   1.31031,   1.30999,   1.30968,   1.30938,   1.30909,
     +   1.30882,   1.30855,   1.30829,   1.30804,   1.30780,   1.30756,
     +   1.30733,   1.30710,   1.30688,   1.30667,   1.30646,   1.30625,
     +   1.30605,   1.30585,   1.30566,   1.30547,   1.30528,   1.30509,
     +   1.30491,   1.30473,   1.30455,   1.30437,   1.30419,   1.30402,
     +   1.30385,   1.30367,   1.30350,   1.30333,   1.30316,   1.30299,
     +   1.30283,   1.30266,   1.30249,   1.30232,   1.30216,   1.30199,
     +   1.30182,   1.30166,   1.30149,   1.30132,   1.30116,   1.30099,
     +   1.30082,   1.30065,   1.30048,   1.30031,   1.30014,   1.29997,
     +   1.29979,   1.29962,   1.29945,   1.29927,   1.29909,   1.29891,
     +   1.29873,   1.29855,   1.29837,   1.29818,   1.29800,   1.29781,
     +   1.29762,   1.29743,   1.29724,   1.29705,   1.29686,   1.29666,
     +   1.29646,   1.29626,   1.29605,   1.29584,   1.29563,   1.29542,
     +   1.29521,   1.29499,   1.29476,   1.29453,   1.29430,   1.29406,
     +   1.29381,   1.29355,   1.29327,   1.29299,   1.29272,   1.29252,
     +   1.29228,   1.29205,   1.29186,   1.29167,   1.29150,   1.29130,
     +   1.29106,   1.29083,   1.29025,   1.28962,   1.28891,   1.28784,
     +   1.28689,   1.28623,   1.28521,   1.28413,   1.28261,   1.28137,
     +   1.28093,   1.28047,   1.28022,   1.27998,   1.27948,   1.27849/
      DATA (TABRE(I),I=229,342)/
     +   1.27774,   1.27691,   1.27610,   1.27535,   1.27471,   1.27404,
     +   1.27329,   1.27240,   1.27139,   1.27029,   1.26901,   1.26736,
     +   1.26591,   1.26441,   1.26284,   1.26036,   1.25860,   1.25815,
     +   1.25768,   1.25675,   1.25579,   1.25383,   1.25179,   1.24967,
     +   1.24745,   1.24512,   1.24266,   1.24004,   1.23725,   1.23270,
     +   1.22583,   1.22198,   1.21548,   1.21184,   1.20790,   1.20507,
     +   1.20209,   1.19566,   1.17411,   1.14734,   1.10766,   1.06739,
     +   1.04762,   1.02650,   1.00357,   0.98197,   0.96503,   0.95962,
     +   0.97269,   0.99172,   1.00668,   1.02186,   1.04270,   1.07597,
     +   1.12954,   1.21267,   1.32509,   1.42599,   1.49656,   1.55095,
     +   1.59988,   1.63631,   1.65024,   1.64278,   1.62691,   1.61284,
     +   1.59245,   1.57329,   1.55770,   1.54129,   1.52654,   1.51139,
     +   1.49725,   1.48453,   1.47209,   1.46125,   1.45132,   1.44215,
     +   1.43366,   1.41553,   1.39417,   1.38732,   1.37735,   1.36448,
     +   1.35414,   1.34456,   1.33882,   1.33807,   1.33847,   1.34053,
     +   1.34287,   1.34418,   1.34634,   1.34422,   1.33453,   1.32897,
     +   1.32333,   1.31800,   1.31432,   1.30623,   1.29722,   1.28898,
     +   1.28730,   1.28603,   1.28509,   1.28535,   1.28813,   1.30156,
     +   1.30901,   1.31720,   1.31893,   1.32039,   1.32201,   1.32239/
      DATA (TABRE(I),I=343,456)/
     +   1.32149,   1.32036,   1.31814,   1.31705,   1.31807,   1.31953,
     +   1.31933,   1.31896,   1.31909,   1.31796,   1.31631,   1.31542,
     +   1.31540,   1.31552,   1.31455,   1.31193,   1.30677,   1.29934,
     +   1.29253,   1.28389,   1.27401,   1.26724,   1.25990,   1.24510,
     +   1.22241,   1.19913,   1.17150,   1.15528,   1.13700,   1.11808,
     +   1.10134,   1.09083,   1.08734,   1.09254,   1.10654,   1.14779,
     +   1.20202,   1.25825,   1.32305,   1.38574,   1.44478,   1.47170,
     +   1.49619,   1.51652,   1.53328,   1.54900,   1.56276,   1.57317,
     +   1.58028,   1.57918,   1.56672,   1.55869,   1.55081,   1.53807,
     +   1.53296,   1.53220,   1.53340,   1.53289,   1.51705,   1.50097,
     +   1.49681,   1.49928,   1.50153,   1.49856,   1.49053,   1.46070,
     +   1.45182,   1.44223,   1.43158,   1.41385,   1.40676,   1.38955,
     +   1.34894,   1.31039,   1.26420,   1.23656,   1.21663,   1.20233,
     +   1.19640,   1.19969,   1.20860,   1.22173,   1.24166,   1.28175,
     +   1.32784,   1.38657,   1.46486,   1.55323,   1.60379,   1.61877,
     +   1.62963,   1.65712,   1.69810,   1.72065,   1.74865,   1.76736,
     +   1.76476,   1.75011,   1.72327,   1.68490,   1.62398,   1.59596,
     +   1.58514,   1.59917,   1.61405,   1.66625,   1.70663,   1.73713,
     +   1.76860,   1.80343,   1.83296,   1.85682,   1.87411,   1.89110/
      DATA (TABRE(I),I=457,468)/
     +   1.89918,   1.90432,   1.90329,   1.88744,   1.87499,   1.86702,
     +   1.85361,   1.84250,   1.83225,   1.81914,   1.82268,   1.82961/
      DATA (TABRET(I,1),I=1,NWLT)/
     +                                    1.82961,   1.83258,   1.83149,
     +   1.82748,   1.82224,   1.81718,   1.81204,   1.80704,   1.80250,
     +   1.79834,   1.79482,   1.79214,   1.78843,   1.78601,   1.78434,
     +   1.78322,   1.78248,   1.78201,   1.78170,   1.78160,   1.78190,
     +   1.78300,   1.78430,   1.78520,   1.78620,   1.78660,   1.78680,
     +   1.78690,   1.78700,   1.78700,   1.78710,   1.78710,   1.78720,
     +   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,
     +   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,
     +   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,
     +   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,
     +   1.78720,   1.78720,   1.78720,   1.78720,   1.78800/
      DATA (TABRET(I,2),I=1,NWLT)/
     +                         1.82961,   1.83258,   1.83149,   1.82748,
     +   1.82224,   1.81718,   1.81204,   1.80704,   1.80250,   1.79834,
     +   1.79482,   1.79214,   1.78843,   1.78601,   1.78434,   1.78322,
     +   1.78248,   1.78201,   1.78170,   1.78160,   1.78190,   1.78300,
     +   1.78430,   1.78520,   1.78610,   1.78630,   1.78640,   1.78650,
     +   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     +   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     +   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     +   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     +   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     +   1.78650,   1.78650,   1.78650,   1.78720/
      DATA(TABRET(I,3),I=1,NWLT)/
     +              1.82961,   1.83258,   1.83149,   1.82748,   1.82224,
     +   1.81718,   1.81204,   1.80704,   1.80250,   1.79834,   1.79482,
     +   1.79214,   1.78843,   1.78601,   1.78434,   1.78322,   1.78248,
     +   1.78201,   1.78160,   1.78140,   1.78160,   1.78220,   1.78310,
     +   1.78380,   1.78390,   1.78400,   1.78400,   1.78400,   1.78400,
     +   1.78400,   1.78390,   1.78380,   1.78370,   1.78370,   1.78370,
     +   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,
     +   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,
     +   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,
     +   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,
     +   1.78370,   1.78400,   1.78450/
      DATA (TABRET(I,4),I=1,NWLT)/
     +   1.82961,   1.83258,   1.83149,   1.82748,   1.82224,   1.81718,
     +   1.81204,   1.80704,   1.80250,   1.79834,   1.79482,   1.79214,
     +   1.78843,   1.78601,   1.78434,   1.78322,   1.78248,   1.78201,
     +   1.78150,   1.78070,   1.78010,   1.77890,   1.77790,   1.77730,
     +   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     +   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     +   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     +   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     +   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     +   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     +   1.77720,   1.77800/
      DATA(TABIM(I),I=1,114)/
     +0.1640E+00,0.1730E+00,0.1830E+00,0.1950E+00,0.2080E+00,0.2230E+00,
     +0.2400E+00,0.2500E+00,0.2590E+00,0.2680E+00,0.2790E+00,0.2970E+00,
     +0.3190E+00,0.3400E+00,0.3660E+00,0.3920E+00,0.4160E+00,0.4400E+00,
     +0.4640E+00,0.4920E+00,0.5170E+00,0.5280E+00,0.5330E+00,0.5340E+00,
     +0.5310E+00,0.5240E+00,0.5100E+00,0.5000E+00,0.4990E+00,0.4680E+00,
     +0.3800E+00,0.3600E+00,0.3390E+00,0.3180E+00,0.2910E+00,0.2510E+00,
     +0.2440E+00,0.2390E+00,0.2390E+00,0.2440E+00,0.2470E+00,0.2240E+00,
     +0.1950E+00,0.1740E+00,0.1720E+00,0.1800E+00,0.1940E+00,0.2130E+00,
     +0.2430E+00,0.2710E+00,0.2890E+00,0.3340E+00,0.3440E+00,0.3820E+00,
     +0.4010E+00,0.4065E+00,0.4050E+00,0.3890E+00,0.3770E+00,0.3450E+00,
     +0.3320E+00,0.3150E+00,0.2980E+00,0.2740E+00,0.2280E+00,0.1980E+00,
     +0.1720E+00,0.1560E+00,0.1100E+00,0.8300E-01,0.5800E-01,0.2200E-01,
     +0.1000E-01,0.3000E-02,0.1000E-02,0.3000E-03,0.1000E-03,0.3000E-04,
     +0.1000E-04,0.3000E-05,0.1000E-05,0.7000E-06,0.4000E-06,0.2000E-06,
     +0.1000E-06,0.6377E-07,0.3750E-07,0.2800E-07,0.2400E-07,0.2200E-07,
     +0.1900E-07,0.1750E-07,0.1640E-07,0.1590E-07,0.1325E-07,0.8623E-08,
     +0.5504E-08,0.3765E-08,0.2710E-08,0.2510E-08,0.2260E-08,0.2080E-08,
     +0.1910E-08,0.1540E-08,0.1530E-08,0.1550E-08,0.1640E-08,0.1780E-08,
     +0.1910E-08,0.2140E-08,0.2260E-08,0.2540E-08,0.2930E-08,0.3110E-08/
      DATA(TABIM(I),I=115,228)/
     +0.3290E-08,0.3520E-08,0.4040E-08,0.4880E-08,0.5730E-08,0.6890E-08,
     +0.8580E-08,0.1040E-07,0.1220E-07,0.1430E-07,0.1660E-07,0.1890E-07,
     +0.2090E-07,0.2400E-07,0.2900E-07,0.3440E-07,0.4030E-07,0.4300E-07,
     +0.4920E-07,0.5870E-07,0.7080E-07,0.8580E-07,0.1020E-06,0.1180E-06,
     +0.1340E-06,0.1400E-06,0.1430E-06,0.1450E-06,0.1510E-06,0.1830E-06,
     +0.2150E-06,0.2650E-06,0.3350E-06,0.3920E-06,0.4200E-06,0.4440E-06,
     +0.4740E-06,0.5110E-06,0.5530E-06,0.6020E-06,0.7550E-06,0.9260E-06,
     +0.1120E-05,0.1330E-05,0.1620E-05,0.2000E-05,0.2250E-05,0.2330E-05,
     +0.2330E-05,0.2170E-05,0.1960E-05,0.1810E-05,0.1740E-05,0.1730E-05,
     +0.1700E-05,0.1760E-05,0.1820E-05,0.2040E-05,0.2250E-05,0.2290E-05,
     +0.3040E-05,0.3840E-05,0.4770E-05,0.5760E-05,0.6710E-05,0.8660E-05,
     +0.1020E-04,0.1130E-04,0.1220E-04,0.1290E-04,0.1320E-04,0.1350E-04,
     +0.1330E-04,0.1320E-04,0.1320E-04,0.1310E-04,0.1320E-04,0.1320E-04,
     +0.1340E-04,0.1390E-04,0.1420E-04,0.1480E-04,0.1580E-04,0.1740E-04,
     +0.1980E-04,0.2500E-04,0.5400E-04,0.1040E-03,0.2030E-03,0.2708E-03,
     +0.3511E-03,0.4299E-03,0.5181E-03,0.5855E-03,0.5899E-03,0.5635E-03,
     +0.5480E-03,0.5266E-03,0.4394E-03,0.3701E-03,0.3372E-03,0.2410E-03,
     +0.1890E-03,0.1660E-03,0.1450E-03,0.1280E-03,0.1030E-03,0.8600E-04,
     +0.8220E-04,0.8030E-04,0.8500E-04,0.9900E-04,0.1500E-03,0.2950E-03/
      DATA(TABIM(I),I=229,342)/
     +0.4687E-03,0.7615E-03,0.1010E-02,0.1313E-02,0.1539E-02,0.1588E-02,
     +0.1540E-02,0.1412E-02,0.1244E-02,0.1068E-02,0.8414E-03,0.5650E-03,
     +0.4320E-03,0.3500E-03,0.2870E-03,0.2210E-03,0.2030E-03,0.2010E-03,
     +0.2030E-03,0.2140E-03,0.2320E-03,0.2890E-03,0.3810E-03,0.4620E-03,
     +0.5480E-03,0.6180E-03,0.6800E-03,0.7300E-03,0.7820E-03,0.8480E-03,
     +0.9250E-03,0.9200E-03,0.8920E-03,0.8700E-03,0.8900E-03,0.9300E-03,
     +0.1010E-02,0.1350E-02,0.3420E-02,0.7920E-02,0.2000E-01,0.3800E-01,
     +0.5200E-01,0.6800E-01,0.9230E-01,0.1270E+00,0.1690E+00,0.2210E+00,
     +0.2760E+00,0.3120E+00,0.3470E+00,0.3880E+00,0.4380E+00,0.4930E+00,
     +0.5540E+00,0.6120E+00,0.6250E+00,0.5930E+00,0.5390E+00,0.4910E+00,
     +0.4380E+00,0.3720E+00,0.3000E+00,0.2380E+00,0.1930E+00,0.1580E+00,
     +0.1210E+00,0.1030E+00,0.8360E-01,0.6680E-01,0.5400E-01,0.4220E-01,
     +0.3420E-01,0.2740E-01,0.2200E-01,0.1860E-01,0.1520E-01,0.1260E-01,
     +0.1060E-01,0.8020E-02,0.6850E-02,0.6600E-02,0.6960E-02,0.9160E-02,
     +0.1110E-01,0.1450E-01,0.2000E-01,0.2300E-01,0.2600E-01,0.2900E-01,
     +0.2930E-01,0.3000E-01,0.2850E-01,0.1730E-01,0.1290E-01,0.1200E-01,
     +0.1250E-01,0.1340E-01,0.1400E-01,0.1750E-01,0.2400E-01,0.3500E-01,
     +0.3800E-01,0.4200E-01,0.4600E-01,0.5200E-01,0.5700E-01,0.6900E-01,
     +0.7000E-01,0.6700E-01,0.6500E-01,0.6400E-01,0.6200E-01,0.5900E-01/
      DATA(TABIM(I),I=343,456)/
     +0.5700E-01,0.5600E-01,0.5500E-01,0.5700E-01,0.5800E-01,0.5700E-01,
     +0.5500E-01,0.5500E-01,0.5400E-01,0.5200E-01,0.5200E-01,0.5200E-01,
     +0.5200E-01,0.5000E-01,0.4700E-01,0.4300E-01,0.3900E-01,0.3700E-01,
     +0.3900E-01,0.4000E-01,0.4200E-01,0.4400E-01,0.4500E-01,0.4600E-01,
     +0.4700E-01,0.5100E-01,0.6500E-01,0.7500E-01,0.8800E-01,0.1080E+00,
     +0.1340E+00,0.1680E+00,0.2040E+00,0.2480E+00,0.2800E+00,0.3410E+00,
     +0.3790E+00,0.4090E+00,0.4220E+00,0.4220E+00,0.4030E+00,0.3890E+00,
     +0.3740E+00,0.3540E+00,0.3350E+00,0.3150E+00,0.2940E+00,0.2710E+00,
     +0.2460E+00,0.1980E+00,0.1640E+00,0.1520E+00,0.1420E+00,0.1280E+00,
     +0.1250E+00,0.1230E+00,0.1160E+00,0.1070E+00,0.7900E-01,0.7200E-01,
     +0.7600E-01,0.7500E-01,0.6700E-01,0.5500E-01,0.4500E-01,0.2900E-01,
     +0.2750E-01,0.2700E-01,0.2730E-01,0.2890E-01,0.3000E-01,0.3400E-01,
     +0.5300E-01,0.7550E-01,0.1060E+00,0.1350E+00,0.1761E+00,0.2229E+00,
     +0.2746E+00,0.3280E+00,0.3906E+00,0.4642E+00,0.5247E+00,0.5731E+00,
     +0.6362E+00,0.6839E+00,0.7091E+00,0.6790E+00,0.6250E+00,0.5654E+00,
     +0.5433E+00,0.5292E+00,0.5070E+00,0.4883E+00,0.4707E+00,0.4203E+00,
     +0.3771E+00,0.3376E+00,0.3056E+00,0.2835E+00,0.3170E+00,0.3517E+00,
     +0.3902E+00,0.4509E+00,0.4671E+00,0.4779E+00,0.4890E+00,0.4899E+00,
     +0.4873E+00,0.4766E+00,0.4508E+00,0.4193E+00,0.3880E+00,0.3433E+00/
      DATA(TABIM(I),I=457,468)/
     +0.3118E+00,0.2935E+00,0.2350E+00,0.1981E+00,0.1865E+00,0.1771E+00,
     +0.1620E+00,0.1490E+00,0.1390E+00,0.1200E+00,0.9620E-01,0.8300E-01/
      DATA(TABIMT(I,1),I=1,NWLT)/
     +                                 0.8300E-01,0.6900E-01,0.5700E-01,
     +0.4560E-01,0.3790E-01,0.3140E-01,0.2620E-01,0.2240E-01,0.1960E-01,
     +0.1760E-01,0.1665E-01,0.1620E-01,0.1550E-01,0.1470E-01,0.1390E-01,
     +0.1320E-01,0.1250E-01,0.1180E-01,0.1060E-01,0.9540E-02,0.8560E-02,
     +0.6210E-02,0.4490E-02,0.3240E-02,0.2340E-02,0.1880E-02,0.1740E-02,
     +0.1500E-02,0.1320E-02,0.1160E-02,0.8800E-03,0.6950E-03,0.4640E-03,
     +0.3400E-03,0.3110E-03,0.2940E-03,0.2790E-03,0.2700E-03,0.2640E-03,
     +0.2580E-03,0.2520E-03,0.2490E-03,0.2540E-03,0.2640E-03,0.2740E-03,
     +0.2890E-03,0.3050E-03,0.3150E-03,0.3460E-03,0.3820E-03,0.4620E-03,
     +0.5000E-03,0.5500E-03,0.5950E-03,0.6470E-03,0.6920E-03,0.7420E-03,
     +0.8200E-03,0.9700E-03,0.1950E-02,0.5780E-02,0.9700E-02/
      DATA(TABIMT(I,2),I=1,NWLT)/
     +                      0.8300E-01,0.6900E-01,0.5700E-01,0.4560E-01,
     +0.3790E-01,0.3140E-01,0.2620E-01,0.2240E-01,0.1960E-01,0.1760E-01,
     +0.1665E-01,0.1600E-01,0.1500E-01,0.1400E-01,0.1310E-01,0.1230E-01,
     +0.1150E-01,0.1080E-01,0.9460E-02,0.8290E-02,0.7270E-02,0.4910E-02,
     +0.3300E-02,0.2220E-02,0.1490E-02,0.1140E-02,0.1060E-02,0.9480E-03,
     +0.8500E-03,0.7660E-03,0.6300E-03,0.5200E-03,0.3840E-03,0.2960E-03,
     +0.2700E-03,0.2520E-03,0.2440E-03,0.2360E-03,0.2300E-03,0.2280E-03,
     +0.2250E-03,0.2200E-03,0.2160E-03,0.2170E-03,0.2200E-03,0.2250E-03,
     +0.2320E-03,0.2390E-03,0.2600E-03,0.2860E-03,0.3560E-03,0.3830E-03,
     +0.4150E-03,0.4450E-03,0.4760E-03,0.5080E-03,0.5400E-03,0.5860E-03,
     +0.6780E-03,0.1280E-02,0.3550E-02,0.5600E-02/
      DATA(TABIMT(I,3),I=1,NWLT)/
     +           0.8300E-01,0.6900E-01,0.5700E-01,0.4560E-01,0.3790E-01,
     +0.3140E-01,0.2620E-01,0.2190E-01,0.1880E-01,0.1660E-01,0.1540E-01,
     +0.1470E-01,0.1350E-01,0.1250E-01,0.1150E-01,0.1060E-01,0.9770E-02,
     +0.9010E-02,0.7660E-02,0.6520E-02,0.5540E-02,0.3420E-02,0.2100E-02,
     +0.1290E-02,0.7930E-03,0.5700E-03,0.5350E-03,0.4820E-03,0.4380E-03,
     +0.4080E-03,0.3500E-03,0.3200E-03,0.2550E-03,0.2120E-03,0.2000E-03,
     +0.1860E-03,0.1750E-03,0.1660E-03,0.1560E-03,0.1490E-03,0.1440E-03,
     +0.1350E-03,0.1210E-03,0.1160E-03,0.1160E-03,0.1170E-03,0.1200E-03,
     +0.1230E-03,0.1320E-03,0.1440E-03,0.1680E-03,0.1800E-03,0.1900E-03,
     +0.2090E-03,0.2160E-03,0.2290E-03,0.2400E-03,0.2600E-03,0.2920E-03,
     +0.6100E-03,0.1020E-02,0.1810E-02/
      DATA(TABIMT(I,4),I=1,NWLT)/
     +0.8300E-01,0.6900E-01,0.5700E-01,0.4450E-01,0.3550E-01,0.2910E-01,
     +0.2440E-01,0.1970E-01,0.1670E-01,0.1400E-01,0.1235E-01,0.1080E-01,
     +0.8900E-02,0.7340E-02,0.6400E-02,0.5600E-02,0.5000E-02,0.4520E-02,
     +0.3680E-02,0.2990E-02,0.2490E-02,0.1550E-02,0.9610E-03,0.5950E-03,
     +0.3690E-03,0.2670E-03,0.2510E-03,0.2290E-03,0.2110E-03,0.1960E-03,
     +0.1730E-03,0.1550E-03,0.1310E-03,0.1130E-03,0.1060E-03,0.9900E-04,
     +0.9300E-04,0.8730E-04,0.8300E-04,0.7870E-04,0.7500E-04,0.6830E-04,
     +0.5600E-04,0.4960E-04,0.4550E-04,0.4210E-04,0.3910E-04,0.3760E-04,
     +0.3400E-04,0.3100E-04,0.2640E-04,0.2510E-04,0.2430E-04,0.2390E-04,
     +0.2370E-04,0.2380E-04,0.2400E-04,0.2460E-04,0.2660E-04,0.4450E-04,
     +0.8700E-04,0.1320E-03/
C
      DATA PI/3.14159265/
C
C     ZERO PARAMETERS
C
      RN=0.0
      CN=0.0
      ABSIND=0.0
      ABSCOF=0.0
C
C     CONVERT WAVELENGTH TO MICRONS
C
      ALAM=XLAM
      IF(IUNIT.EQ.1)ALAM=1000*ALAM
      IF(IUNIT.EQ.2)ALAM=10000*ALAM
      IF(IUNIT.EQ.3)ALAM=10000*(1.0/ALAM)
      IF(ALAM.LT.WLMIN.OR.ALAM.GT.WLMAX)RETURN
      IF(ALAM.GT.CUTICE)GO TO 10
C
C     REGION FROM 0.045 MICRONS TO 167.0 MICRONS - NO TEMPERATURE DEPEND
C
      DO 1 I=2,NWL
      IF(ALAM.LT.WL(I)) GO TO 2
 1    CONTINUE
 2    X1=ALOG(WL(I-1))
      X2=ALOG(WL(I))
      Y1=TABRE(I-1)
      Y2=TABRE(I)
      X=ALOG(ALAM)
      Y=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
      RN=Y
      Y1=ALOG(ABS(TABIM(I-1)))
      Y2=ALOG(ABS(TABIM(I)))
      Y=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
      CN=EXP(Y)
      GO TO 20
C
C     REGION FROM 167.0 MICRONS TO 8.6 METERS - TEMPERATURE DEPENDENCE
C
 10   TK=T
      IF(TK.GT.TEMREF(1))TK=TEMREF(1)
      IF(TK.LT.TEMREF(4))TK=TEMREF(4)
      DO 11 I=2,4
      IF(TK.GE.TEMREF(I)) GO TO 12
 11   CONTINUE
 12   LT1=I
      LT2=I-1
      DO 13 I=2,NWLT
      IF(ALAM.LE.WLT(I)) GO TO 14
 13   CONTINUE
 14   X1=ALOG(WLT(I-1))
      X2=ALOG(WLT(I))
      Y1=TABRET(I-1,LT1)
      Y2=TABRET(I,LT1)
      X=ALOG(ALAM)
      YLO=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
      Y1=TABRET(I-1,LT2)
      Y2=TABRET(I,LT2)
      YHI=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
      T1=TEMREF(LT1)
      T2=TEMREF(LT2)
      Y=((TK-T1)*(YHI-YLO)/(T2-T1))+YLO
      RN=Y
      Y1=ALOG(ABS(TABIMT(I-1,LT1)))
      Y2=ALOG(ABS(TABIMT(I,LT1)))
      YLO=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
      Y1=ALOG(ABS(TABIMT(I-1,LT2)))
      Y2=ALOG(ABS(TABIMT(I,LT2)))
      YHI=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
      Y=((TK-T1)*(YHI-YLO)/(T2-T1))+YLO
      CN=EXP(Y)
C
C     ABSORPTIVE QUANITIES
C
 20   ABSIND=CN/RN
      ABSCOF=4.0*PI*CN/ALAM
      RETURN
      END

      SUBROUTINE TABLE2(X,Y,Z)
      IMPLICIT NONE

C Parameters:

      INTEGER NTABMX
      PARAMETER(NTABMX=10000)

C Arguments:

      REAL X,Y,Z

C Common variables:

      INTEGER NTAB
      REAL XT,YT,ZT
      COMMON/TABCOM2/XT(1:NTABMX),YT(1:NTABMX),ZT(1:NTABMX),NTAB

c Local variables:

      INTEGER I1,I2,I3,I4,N2,NOUT
      REAL X1,X2,X3,X4,Y1,Y2,Y3,Y4

c***********************************************************************
c Subroutine TABLE2

c Given:
c       X
c       XT(1-NTAB),YT(1-NTAB),ZT(1-NTAB),NTAB (in COMMON/TABCOM2/)
c Returns:
c       Y = y(X), Z=z(X) where y(X) is local parabolic fit to {XT,YT}
c                        z(X) is local parabolic fit to {XT,ZT}
c
c If X.LE.XT(2) or X.GE.XT(NTAB-1) then use parabola passing
c through first 3 or last 3 points.
c Otherwise, use parabola passing through points just below and
c above X, and minimizing sum of squared deviations to next point
c on either side.
c
c B.T.Draine, I.A.S., March 1980
c History:
c 08.03.01 (BTD): first written at IAS
c 91.07.21 (BTD): explicit declaration of all variables
c 99.02.17 (BTD): add SAVE statements to make f77 compliant and
c                 g77 compatible
c 05.01.21 (BTD): major rewrite
c 07.11.16 (BTD): add more output when detect outside table limits
c 07.11.16 (BTD): corrected typo: was not calculating Z 
c                 when near table boundary (inadvertently had Y in place
c                 of Z in call to PARAB3 )
c 11.06.14 (BTD): change
c                 IF((X-XT(2))*(XT(NTAB-1)-X).GE.0.D0)THEN
c                 to
c                 IF((X-XT(2))*(XT(NTAB-1)-X).GT.0.D0)THEN
c 11.07.11 (BTD): change
c                 IF(X.LT.XT(2))THEN
c                 to
c                 IF(X.LE.XT(2))THEN
c                 change
c                 IF(X.GT.XT(2))THEN
c                 to
c                 IF(X.GE.XT(2))THEN
c end history
c***********************************************************************

      DATA I2/2/,N2/1/,NOUT/0/
      SAVE I2,N2,NOUT

c sanity check:

      IF(NTAB.LT.4)THEN
         IF(NTAB.EQ.3)THEN
            X1=XT(1)
            X2=XT(2)
            X3=XT(3)
            Y1=YT(1)
            Y2=YT(2)
            Y3=YT(3)
c*** diagnostic
c            write(0,*)'table2 ckpt 2, x=',x
c***
            CALL PARAB3(X,Y,X1,X2,X3,Y1,Y2,Y3)
            Y1=ZT(1)
            Y2=ZT(2)
            Y3=ZT(3)
            CALL PARAB3(X,Z,X1,X2,X3,Y1,Y2,Y3)
            IF((X-X1)*(X2-X).LT.0.)THEN
               NOUT=NOUT+1
               WRITE(0,6999)X,NOUT
               N2=2*N2
            ENDIF
            RETURN
         ELSE
            WRITE(0,*)' FATAL ERROR: attempt to use subroutine',
     &                ' TABLE2 for array with < 3 elements'
            STOP
         ENDIF
      ENDIF

c we can assume at least 4 elements in array
c is X between XT(2) and XT(NTAB-1) ?

      IF((X-XT(2))*(XT(NTAB-1)-X).GT.0.D0)THEN

c not near table boundaries: use 4 point interpolation
c determine I1,I2,I3,I4

         CALL BRACKET(X,XT,NTAB,I2)

         I1=I2-1
         I3=I2+1
         I4=I2+2
c*** sanity check
         if(i1.lt.1)then
            write(0,*)'table2 ckpt 3: i1=',i1,' < 1'
            write(0,*)'  xt(2)   =',xt(2)
            write(0,*)'    x     =',x
            write(0,*)'xt(ntab-1)=',xt(ntab-1)
            write(0,*)'fatal error: stop'
            stop
         elseif(i4.gt.ntab)then
            write(0,*)'table2 ckpt 4: i4=',i4,' > ntab'
            write(0,*)'fatal error: stop'
            stop
         endif
c***
c*** diagnostic
c         write(0,*)'table2 ckpt 5: x=',x
c         write(0,*)'   i1,i2,i3,i4=',i1,i2,i3,i4
c         write(0,*)'   xt(i1)=',xt(i1)
c         write(0,*)'   xt(i2)=',xt(i2)
c         write(0,*)'   xt(i3)=',xt(i3)
c         write(0,*)'   xt(i4)=',xt(i4)
c***
         X1=XT(I1)
         X2=XT(I2)
         X3=XT(I3)
         X4=XT(I4)
         Y1=YT(I1)
         Y2=YT(I2)
         Y3=YT(I3)
         Y4=YT(I4)
c*** diagnostic
c         write(0,*)' yt(i1)=',yt(i1)
c         write(0,*)' yt(i2)=',yt(i2)
c         write(0,*)' yt(i3)=',yt(i3)
c         write(0,*)' yt(i4)=',yt(i4)
c***
         CALL PARAB4(X,Y,X1,X2,X3,X4,Y1,Y2,Y3,Y4)
c*** diagnostic
c         write(0,*)'    y  =',y
c         write(0,*)' zt(i1)=',zt(i1)
c         write(0,*)' zt(i2)=',zt(i2)
c         write(0,*)' zt(i3)=',zt(i3)
c         write(0,*)' zt(i4)=',zt(i4)
c***
         Y1=ZT(I1)
         Y2=ZT(I2)
         Y3=ZT(I3)
         Y4=ZT(I4)
         CALL PARAB4(X,Z,X1,X2,X3,X4,Y1,Y2,Y3,Y4)
c*** diagnostic
c         write(0,*)'    z  =',z
c***
      ELSE

c near table boundary: fix I1,I2,I3

         IF(XT(NTAB).GT.XT(1))THEN

c*** diagnostic
c            write(0,*)'table2 ckpt 6, x=',x
c***
c x monotonically increasing

c 110711 BTD change
c            IF(X.LT.XT(2))THEN
            IF(X.LE.XT(2))THEN
               I1=1
               I2=2
               I3=3
               X3=XT(3)
            ELSE
               I1=NTAB-2
               I2=NTAB-1
               I3=NTAB
            ENDIF
         ELSE

c x monotonically decreasing

c 110711 BTD change
c            IF(X.GT.XT(2))THEN
            IF(X.GE.XT(2))THEN
               I1=1
               I2=2
               I3=3
            ELSE
               I1=NTAB-2
               I2=NTAB-1
               I3=NTAB
            ENDIF
         ENDIF

c check to see if outside table limits. if so, issue warning

         IF((X-XT(1))*(XT(NTAB)-X).LT.0.D0)THEN
            NOUT=NOUT+1
            IF(NOUT/N2.EQ.1)THEN
               WRITE(0,6999)X,NOUT
               WRITE(0,7999)xt(1),xt(ntab),ntab
               N2=2*N2
            ENDIF
         ENDIF
         X1=XT(I1)
         X2=XT(I2)
         X3=XT(I3)
         Y1=YT(I1)
         Y2=YT(I2)
         Y3=YT(I3)
c*** diagnostic
c         write(0,*)'table2 ckpt 7, x=',x
c         write(0,*)'   x1=',x1,' y1=',y1
c         write(0,*)'   x2=',x2,' y2=',y2
c         write(0,*)'   x3=',x3,' y3=',y3
c***
         CALL PARAB3(X,Y,X1,X2,X3,Y1,Y2,Y3)
c*** diagnostic
c         write(0,*)'   y=',y
c***
         Y1=ZT(I1)
         Y2=ZT(I2)
         Y3=ZT(I3)
c*** diagnostic
c         write(0,*)'table2 ckpt 8, x=',x
c         write(0,*)'   x1=',x1,' y1=',y1
c         write(0,*)'   x2=',x2,' y2=',y2
c         write(0,*)'   x3=',x3,' y3=',y3
c***
         CALL PARAB3(X,Z,X1,X2,X3,Y1,Y2,Y3)
c*** diagnostic
c         write(0,*)'   z=',z
c***
      ENDIF
      RETURN
 6999 FORMAT(' WARNING: outside table limits for X=',1PE12.5,
     &       ' NOUT=',I5)
 7999 FORMAT(' XT(1)=',1pe12.5,' xt(ntab)=',e12.5,' ntab=',I6)
      END

      SUBROUTINE BRACKET(X,XT,NTAB,I)
      IMPLICIT NONE

c arguments:

      INTEGER I,NTAB
      REAL X
      REAL XT(1:NTAB)

c local variables:

      INTEGER I1,I2,INC,ITRY

c=======================================================================
c subroutine BRACKET
c given:
c     X    = real variable
c     XT(1-NTAB)=array, must be monotonically-increasing
c     NTAB = number of elements in array
c     I    = index that may be near desired I

c returns
c
c     I = index such that X lies between XT(I) and XT(I+1)
c         if X < XT(1) : Issue error message and stop
c         if XT(1).LE.X.LE.XT(2), then set I=2
c         if XT(NTAB-1).LT.X.LE.XT(NTAB), then set I=NTAB-1
c         if X > XT(NTAB): Issue error message and stop
c
c B.T. Draine, Princeton University, 2005.01.21
c history
c 05.01.21 (BTD): first written, to replace old code in TABLE
c 05.06.23 (BTD): rewrote to correct error
C end history
c=======================================================================
c sanity checks

      IF(XT(NTAB).LE.XT(1))THEN
         WRITE(0,*)'Fatal error in BRACKET: array not ',
     &             'monotonically-increasing'
         STOP
      ENDIF
      IF((X-XT(1))*(XT(NTAB)-X).LT.0.)THEN
         WRITE(0,*)'Fatal error: BRACKET called with X=',X,
     &             'outside XT(1),XT(NTAB)=',XT(1),XT(NTAB)
         STOP
      ENDIF

c special treatment for boundaries:

      IF(X.LT.XT(2))THEN
         I=2
         RETURN
      ELSEIF(X.GT.XT(NTAB-1))THEN
         I=NTAB-1
         RETURN
      ENDIF

c general treatment:
c find I1 and I2 bracketing x

      INC=1

c if input I is in allowed range, use it as starting point

      IF(I.GT.0..AND.I.LE.NTAB)THEN
         I1=I
      ELSE
         I1=(1+NTAB)/2
      ENDIF
      I2=I

 0100 I1=MAX(I1-INC,1)
      I2=MIN(I2+INC,NTAB)
      IF((X-XT(I1))*(XT(I2)-X).LT.0.)THEN
         INC=INC*2
         GOTO 0100
      ENDIF

c have coarsely bracketed I
c use binary chop to find I
      
 1000 IF(I2-I1.GT.1)THEN
         ITRY=I1+(I2-I1)/2
         IF((X-XT(I1))*(XT(ITRY)-X).GE.0)THEN

c in interval [i1,itry]

            I2=ITRY
            GOTO 1000
         ELSE

c in interval (itry,i2]

            I1=ITRY
            GOTO 1000
         ENDIF
      ENDIF
      I=I1
c*** diagnostic
      if(i.lt.1)then
         write(0,*)'table 2 ckpt 9, i=',i
      endif
c***
      RETURN
      END
      
      SUBROUTINE PARAB3(X,Y,X1,X2,X3,Y1,Y2,Y3)
      IMPLICIT NONE

C Arguments:

      REAL X,X1,X2,X3,Y,Y1,Y2,Y3

C Local variables:

      REAL A,B

C=======================================================================
C Subroutine PARAB3 
C Given:
C    X
C    X1,Y1
C    X2,Y2
C    X3,Y3
C Returns
C    Y = f(X), where f(X) is parabola
C        constrained to fit (X1,Y1),(X2,Y2),(X3,Y3) exactly

C B.T.Draine, I.A.S., March 1980

C history:
C 91.07.21 (BTD): explicit declaration of all variables
C 97.02.17 (BTD): add SAVE statement to make f77 compliant and
C                 g77 compatible
C 05.01.21 (BTD): revised to no longer attempt to use previously-
C                 computed values of A and B.
C end history
C***********************************************************************
      A=(Y3-Y2-(X3-X2)*(Y1-Y2)/(X1-X2))/((X3-X2)*(X3-X1))
      B=(Y1-Y2)/(X1-X2)-A*(X1-X2)
      Y=(A*(X-X2)+B)*(X-X2)+Y2
      RETURN
      END

      SUBROUTINE PARAB4(X,Y,X1,X2,X3,X4,Y1,Y2,Y3,Y4)
      IMPLICIT NONE

C Arguments:

      REAL X,X1,X2,X3,X4,Y,Y1,Y2,Y3,Y4

C Local variables:

      REAL A,B

C=======================================================================
C Subroutine PARAB4 
C Given:
C    X
C    X1,Y1
C    X2,Y2
C    X3,Y3
C    X4,Y4
C Returns:
C    Y = f(X), where f(X) is parabola constrained to pass through 
C        (X2,Y2) and (X3,Y3), and to minimize sum of squared deviation 
C        from (X1,Y1) and (X4,Y4).
C
C    It is assumed that X1.LT.X2.LE.X.LT.X3.LT.X4
C    Parabola is given by Y=A*(X-X2)**2+B*(X-X2)+Y2

C B.T. Draine, I.A.S., March 1980
C History:
C 91.07.21 (BTD): explicit declaration of all variables.
C 97.02.17 (BTD): add SAVE statement to make f77 compliant and
C                 g77 compatible
C 05.01.21 (BTD): revised to no longer attempt to use previously-
C                 computed values of A and B
C end history
C***********************************************************************
      A=((X1-X2)*(Y3-Y2)/(X3-X2)+Y2-Y1)*(X1-X2)*(X1-X3)
     & +((X4-X2)*(Y3-Y2)/(X3-X2)+Y2-Y4)*(X4-X2)*(X4-X3)
      A=-A/(((X1-X2)*(X1-X3))**2+((X4-X2)*(X4-X3))**2)
      B=(Y3-Y2)/(X3-X2)-A*(X3-X2)
      Y=(A*(X-X2)+B)*(X-X2)+Y2
      RETURN
      END
*DECK XERCNT
      SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)
C***BEGIN PROLOGUE  XERCNT
C***SUBSIDIARY
C***PURPOSE  Allow user control over handling of errors.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERCNT-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        Allows user control over handling of individual errors.
C        Just after each message is recorded, but before it is
C        processed any further (i.e., before it is printed or
C        a decision to abort is made), a call is made to XERCNT.
C        If the user has provided his own version of XERCNT, he
C        can then override the value of KONTROL used in processing
C        this message by redefining its value.
C        KONTRL may be set to any value from -2 to 2.
C        The meanings for KONTRL are the same as in XSETF, except
C        that the value of KONTRL changes only for this message.
C        If KONTRL is set to a value outside the range from -2 to 2,
C        it will be moved back into that range.
C
C     Description of Parameters
C
C      --Input--
C        LIBRAR - the library that the routine is in.
C        SUBROU - the subroutine that XERMSG is being called from
C        MESSG  - the first 20 characters of the error message.
C        NERR   - same as in the call to XERMSG.
C        LEVEL  - same as in the call to XERMSG.
C        KONTRL - the current value of the control flag as set
C                 by a call to XSETF.
C
C      --Output--
C        KONTRL - the new value of KONTRL.  If KONTRL is not
C                 defined, it will remain at its original value.
C                 This changed value of control affects only
C                 the current occurrence of the current message.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900206  Routine changed from user-callable to subsidiary.  (WRB)
C   900510  Changed calling sequence to include LIBRARY and SUBROUTINE
C           names, changed routine name from XERCTL to XERCNT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERCNT
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
C***FIRST EXECUTABLE STATEMENT  XERCNT
      RETURN
      END
*DECK XERHLT
      SUBROUTINE XERHLT (MESSG)
C***BEGIN PROLOGUE  XERHLT
C***SUBSIDIARY
C***PURPOSE  Abort program execution and print error message.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERHLT-A)
C***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        ***Note*** machine dependent routine
C        XERHLT aborts the execution of the program.
C        The error message causing the abort is given in the calling
C        sequence, in case one needs it for printing on a dayfile,
C        for example.
C
C     Description of Parameters
C        MESSG is as in XERMSG.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900206  Routine changed from user-callable to subsidiary.  (WRB)
C   900510  Changed calling sequence to delete length of character
C           and changed routine name from XERABT to XERHLT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERHLT
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERHLT
      STOP
      END
*DECK XERMSG
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
C***BEGIN PROLOGUE  XERMSG
C***PURPOSE  Process error messages for SLATEC and other libraries.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERMSG-A)
C***KEYWORDS  ERROR MESSAGE, XERROR
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
C***DESCRIPTION
C
C   XERMSG processes a diagnostic message in a manner determined by the
C   value of LEVEL and the current value of the library error control
C   flag, KONTRL.  See subroutine XSETF for details.
C
C    LIBRAR   A character constant (or character variable) with the name
C             of the library.  This will be 'SLATEC' for the SLATEC
C             Common Math Library.  The error handling package is
C             general enough to be used by many libraries
C             simultaneously, so it is desirable for the routine that
C             detects and reports an error to identify the library name
C             as well as the routine name.
C
C    SUBROU   A character constant (or character variable) with the name
C             of the routine that detected the error.  Usually it is the
C             name of the routine that is calling XERMSG.  There are
C             some instances where a user callable library routine calls
C             lower level subsidiary routines where the error is
C             detected.  In such cases it may be more informative to
C             supply the name of the routine the user called rather than
C             the name of the subsidiary routine that detected the
C             error.
C
C    MESSG    A character constant (or character variable) with the text
C             of the error or warning message.  In the example below,
C             the message is a character constant that contains a
C             generic message.
C
C                   CALL XERMSG ('SLATEC', 'MMPY',
C                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
C                  *3, 1)
C
C             It is possible (and is sometimes desirable) to generate a
C             specific message--e.g., one that contains actual numeric
C             values.  Specific numeric values can be converted into
C             character strings using formatted WRITE statements into
C             character variables.  This is called standard Fortran
C             internal file I/O and is exemplified in the first three
C             lines of the following example.  You can also catenate
C             substrings of characters to construct the error message.
C             Here is an example showing the use of both writing to
C             an internal file and catenating character strings.
C
C                   CHARACTER*5 CHARN, CHARL
C                   WRITE (CHARN,10) N
C                   WRITE (CHARL,10) LDA
C                10 FORMAT(I5)
C                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
C                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
C                  *   CHARL, 3, 1)
C
C             There are two subtleties worth mentioning.  One is that
C             the // for character catenation is used to construct the
C             error message so that no single character constant is
C             continued to the next line.  This avoids confusion as to
C             whether there are trailing blanks at the end of the line.
C             The second is that by catenating the parts of the message
C             as an actual argument rather than encoding the entire
C             message into one large character variable, we avoid
C             having to know how long the message will be in order to
C             declare an adequate length for that large character
C             variable.  XERMSG calls XERPRN to print the message using
C             multiple lines if necessary.  If the message is very long,
C             XERPRN will break it into pieces of 72 characters (as
C             requested by XERMSG) for printing on multiple lines.
C             Also, XERMSG asks XERPRN to prefix each line with ' *  '
C             so that the total line length could be 76 characters.
C             Note also that XERPRN scans the error message backwards
C             to ignore trailing blanks.  Another feature is that
C             the substring '$$' is treated as a new line sentinel
C             by XERPRN.  If you want to construct a multiline
C             message without having to count out multiples of 72
C             characters, just use '$$' as a separator.  '$$'
C             obviously must occur within 72 characters of the
C             start of each line to have its intended effect since
C             XERPRN is asked to wrap around at 72 characters in
C             addition to looking for '$$'.
C
C    NERR     An integer value that is chosen by the library routine's
C             author.  It must be in the range -99 to 999 (three
C             printable digits).  Each distinct error should have its
C             own error number.  These error numbers should be described
C             in the machine readable documentation for the routine.
C             The error numbers need be unique only within each routine,
C             so it is reasonable for each routine to start enumerating
C             errors from 1 and proceeding to the next integer.
C
C    LEVEL    An integer value in the range 0 to 2 that indicates the
C             level (severity) of the error.  Their meanings are
C
C            -1  A warning message.  This is used if it is not clear
C                that there really is an error, but the user's attention
C                may be needed.  An attempt is made to only print this
C                message once.
C
C             0  A warning message.  This is used if it is not clear
C                that there really is an error, but the user's attention
C                may be needed.
C
C             1  A recoverable error.  This is used even if the error is
C                so serious that the routine cannot return any useful
C                answer.  If the user has told the error package to
C                return after recoverable errors, then XERMSG will
C                return to the Library routine which can then return to
C                the user's routine.  The user may also permit the error
C                package to terminate the program upon encountering a
C                recoverable error.
C
C             2  A fatal error.  XERMSG will not return to its caller
C                after it receives a fatal error.  This level should
C                hardly ever be used; it is much better to allow the
C                user a chance to recover.  An example of one of the few
C                cases in which it is permissible to declare a level 2
C                error is a reverse communication Library routine that
C                is likely to be called repeatedly until it integrates
C                across some interval.  If there is a serious error in
C                the input such that another step cannot be taken and
C                the Library routine is called again without the input
C                error having been corrected by the caller, the Library
C                routine will probably be called forever with improper
C                input.  In this case, it is reasonable to declare the
C                error to be fatal.
C
C    Each of the arguments to XERMSG is input; none will be modified by
C    XERMSG.  A routine may make multiple calls to XERMSG with warning
C    level messages; however, after a call to XERMSG with a recoverable
C    error, the routine should return to the user.  Do not try to call
C    XERMSG with a second recoverable error after the first recoverable
C    error because the error package saves the error number.  The user
C    can retrieve this error number by calling another entry point in
C    the error handling package and then clear the error number when
C    recovering from the error.  Calling XERMSG in succession causes the
C    old error number to be overwritten by the latest error number.
C    This is considered harmless for error numbers associated with
C    warning messages but must not be done for error numbers of serious
C    errors.  After a call to XERMSG with a recoverable error, the user
C    must be given a chance to call NUMXER or XERCLR to retrieve or
C    clear the error number.
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
C***REVISION HISTORY  (YYMMDD)
C   880101  DATE WRITTEN
C   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
C           THERE ARE TWO BASIC CHANGES.
C           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
C               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
C               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
C               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
C               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
C               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
C               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
C               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
C           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
C               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
C               OF LOWER CASE.
C   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
C           THE PRINCIPAL CHANGES ARE
C           1.  CLARIFY COMMENTS IN THE PROLOGUES
C           2.  RENAME XRPRNT TO XERPRN
C           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
C               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
C               CHARACTER FOR NEW RECORDS.
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
C           CLEAN UP THE CODING.
C   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
C           PREFIX.
C   891013  REVISED TO CORRECT COMMENTS.
C   891214  Prologue converted to Version 4.0 format.  (WRB)
C   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
C           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
C           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
C           XERCTL to XERCNT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERMSG
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8 XLIBR, XSUBR
      CHARACTER*72  TEMP
      CHARACTER*20  LFIRST
C***FIRST EXECUTABLE STATEMENT  XERMSG
      LKNTRL = J4SAVE (2, 0, .FALSE.)
      MAXMES = J4SAVE (4, 0, .FALSE.)
C
C       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
C       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
C          SHOULD BE PRINTED.
C
C       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
C          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
C          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
C
      IF (NERR.LT.-9999999 .OR. NERR.GT.99999999 .OR. NERR.EQ.0 .OR.
     *   LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN
         CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' //
     *      'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '//
     *      'JOB ABORT DUE TO FATAL ERROR.', 72)
         CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
         CALL XERHLT (' ***XERMSG -- INVALID INPUT')
         RETURN
      ENDIF
C
C       RECORD THE MESSAGE.
C
      I = J4SAVE (1, NERR, .TRUE.)
      CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
C
C       HANDLE PRINT-ONCE WARNING MESSAGES.
C
      IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN
C
C       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
C
      XLIBR  = LIBRAR
      XSUBR  = SUBROU
      LFIRST = MESSG
      LERR   = NERR
      LLEVEL = LEVEL
      CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
C
      LKNTRL = MAX(-2, MIN(2,LKNTRL))
      MKNTRL = ABS(LKNTRL)
C
C       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
C       ZERO AND THE ERROR IS NOT FATAL.
C
      IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30
      IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30
      IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30
      IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30
C
C       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
C       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
C       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
C       IS NOT ZERO.
C
      IF (LKNTRL .NE. 0) THEN
         TEMP(1:21) = 'MESSAGE FROM ROUTINE '
         I = MIN(LEN(SUBROU), 16)
         TEMP(22:21+I) = SUBROU(1:I)
         TEMP(22+I:33+I) = ' IN LIBRARY '
         LTEMP = 33 + I
         I = MIN(LEN(LIBRAR), 16)
         TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
         TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
         LTEMP = LTEMP + I + 1
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
C       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
C       FROM EACH OF THE FOLLOWING THREE OPTIONS.
C       1.  LEVEL OF THE MESSAGE
C              'INFORMATIVE MESSAGE'
C              'POTENTIALLY RECOVERABLE ERROR'
C              'FATAL ERROR'
C       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
C              'PROG CONTINUES'
C              'PROG ABORTED'
C       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
C           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
C           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
C              'TRACEBACK REQUESTED'
C              'TRACEBACK NOT REQUESTED'
C       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
C       EXCEED 74 CHARACTERS.
C       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
C
      IF (LKNTRL .GT. 0) THEN
C
C       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
C
         IF (LEVEL .LE. 0) THEN
            TEMP(1:20) = 'INFORMATIVE MESSAGE,'
            LTEMP = 20
         ELSEIF (LEVEL .EQ. 1) THEN
            TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
            LTEMP = 30
         ELSE
            TEMP(1:12) = 'FATAL ERROR,'
            LTEMP = 12
         ENDIF
C
C       THEN WHETHER THE PROGRAM WILL CONTINUE.
C
         IF ((MKNTRL.EQ.2 .AND. LEVEL.GE.1) .OR.
     *       (MKNTRL.EQ.1 .AND. LEVEL.EQ.2)) THEN
            TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
            LTEMP = LTEMP + 14
         ELSE
            TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
            LTEMP = LTEMP + 16
         ENDIF
C
C       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
C
         IF (LKNTRL .GT. 0) THEN
            TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
            LTEMP = LTEMP + 20
         ELSE
            TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
            LTEMP = LTEMP + 24
         ENDIF
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       NOW SEND OUT THE MESSAGE.
C
      CALL XERPRN (' *  ', -1, MESSG, 72)
C
C       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
C          TRACEBACK.
C
      IF (LKNTRL .GT. 0) THEN
         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
         DO 10 I=16,22
            IF (TEMP(I:I) .NE. ' ') GO TO 20
   10    CONTINUE
C
   20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
         CALL FDUMP
      ENDIF
C
C       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
C
      IF (LKNTRL .NE. 0) THEN
         CALL XERPRN (' *  ', -1, ' ', 72)
         CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)
         CALL XERPRN ('    ',  0, ' ', 72)
      ENDIF
C
C       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
C       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
C
   30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN
C
C       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
C       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
C       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
C
      IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN
         IF (LEVEL .EQ. 1) THEN
            CALL XERPRN
     *         (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
         ELSE
            CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
         ENDIF
         CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
         CALL XERHLT (' ')
      ELSE
         CALL XERHLT (MESSG)
      ENDIF
      RETURN
      END
*DECK XERPRN
      SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP)
C***BEGIN PROLOGUE  XERPRN
C***SUBSIDIARY
C***PURPOSE  Print error messages processed by XERMSG.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERPRN-A)
C***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
C***DESCRIPTION
C
C This routine sends one or more lines to each of the (up to five)
C logical units to which error messages are to be sent.  This routine
C is called several times by XERMSG, sometimes with a single line to
C print and sometimes with a (potentially very long) message that may
C wrap around into multiple lines.
C
C PREFIX  Input argument of type CHARACTER.  This argument contains
C         characters to be put at the beginning of each line before
C         the body of the message.  No more than 16 characters of
C         PREFIX will be used.
C
C NPREF   Input argument of type INTEGER.  This argument is the number
C         of characters to use from PREFIX.  If it is negative, the
C         intrinsic function LEN is used to determine its length.  If
C         it is zero, PREFIX is not used.  If it exceeds 16 or if
C         LEN(PREFIX) exceeds 16, only the first 16 characters will be
C         used.  If NPREF is positive and the length of PREFIX is less
C         than NPREF, a copy of PREFIX extended with blanks to length
C         NPREF will be used.
C
C MESSG   Input argument of type CHARACTER.  This is the text of a
C         message to be printed.  If it is a long message, it will be
C         broken into pieces for printing on multiple lines.  Each line
C         will start with the appropriate prefix and be followed by a
C         piece of the message.  NWRAP is the number of characters per
C         piece; that is, after each NWRAP characters, we break and
C         start a new line.  In addition the characters '$$' embedded
C         in MESSG are a sentinel for a new line.  The counting of
C         characters up to NWRAP starts over for each new line.  The
C         value of NWRAP typically used by XERMSG is 72 since many
C         older error messages in the SLATEC Library are laid out to
C         rely on wrap-around every 72 characters.
C
C NWRAP   Input argument of type INTEGER.  This gives the maximum size
C         piece into which to break MESSG for printing on multiple
C         lines.  An embedded '$$' ends a line, and the count restarts
C         at the following character.  If a line break does not occur
C         on a blank (it would split a word) that word is moved to the
C         next line.  Values of NWRAP less than 16 will be treated as
C         16.  Values of NWRAP greater than 132 will be treated as 132.
C         The actual line length will be NPREF + NWRAP after NPREF has
C         been adjusted to fall between 0 and 16 and NWRAP has been
C         adjusted to fall between 16 and 132.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  I1MACH, XGETUA
C***REVISION HISTORY  (YYMMDD)
C   880621  DATE WRITTEN
C   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
C           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
C           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
C           SLASH CHARACTER IN FORMAT STATEMENTS.
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
C           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
C           LINES TO BE PRINTED.
C   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
C           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
C   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
C   891214  Prologue converted to Version 4.0 format.  (WRB)
C   900510  Added code to break messages between words.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERPRN
      CHARACTER*(*) PREFIX, MESSG
      INTEGER NPREF, NWRAP
      CHARACTER*148 CBUFF
      INTEGER IU(5), NUNIT
      CHARACTER*2 NEWLIN
      PARAMETER (NEWLIN = '$$')
C***FIRST EXECUTABLE STATEMENT  XERPRN
      CALL XGETUA(IU,NUNIT)
C
C       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
C       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
C       ERROR MESSAGE UNIT.
C
      N = I1MACH(4)
      DO 10 I=1,NUNIT
         IF (IU(I) .EQ. 0) IU(I) = N
   10 CONTINUE
C
C       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
C       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
C       THE REST OF THIS ROUTINE.
C
      IF ( NPREF .LT. 0 ) THEN
         LPREF = LEN(PREFIX)
      ELSE
         LPREF = NPREF
      ENDIF
      LPREF = MIN(16, LPREF)
      IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX
C
C       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
C       TIME FROM MESSG TO PRINT ON ONE LINE.
C
      LWRAP = MAX(16, MIN(132, NWRAP))
C
C       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
C
      LENMSG = LEN(MESSG)
      N = LENMSG
      DO 20 I=1,N
         IF (MESSG(LENMSG:LENMSG) .NE. ' ') GO TO 30
         LENMSG = LENMSG - 1
   20 CONTINUE
   30 CONTINUE
C
C       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
C
      IF (LENMSG .EQ. 0) THEN
         CBUFF(LPREF+1:LPREF+1) = ' '
         DO 40 I=1,NUNIT
            WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
   40    CONTINUE
         RETURN
      ENDIF
C
C       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
C       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
C       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
C       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
C
C       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
C       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
C       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
C       OF THE SECOND ARGUMENT.
C
C       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
C       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
C       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
C       POSITION NEXTC.
C
C       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
C                       REMAINDER OF THE CHARACTER STRING.  LPIECE
C                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
C                       WHICHEVER IS LESS.
C
C       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
C                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
C                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
C                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
C                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
C                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
C                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
C                       SHOULD BE INCREMENTED BY 2.
C
C       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.
C
C       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1
C                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
C                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ.
C                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
C                       AT THE END OF A LINE.
C
      NEXTC = 1
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
      IF (LPIECE .EQ. 0) THEN
C
C       THERE WAS NO NEW LINE SENTINEL FOUND.
C
         IDELTA = 0
         LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
         IF (LPIECE .LT. LENMSG+1-NEXTC) THEN
            DO 52 I=LPIECE+1,2,-1
               IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
                  LPIECE = I-1
                  IDELTA = 1
                  GOTO 54
               ENDIF
   52       CONTINUE
         ENDIF
   54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSEIF (LPIECE .EQ. 1) THEN
C
C       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
C       DON'T PRINT A BLANK LINE.
C
         NEXTC = NEXTC + 2
         GO TO 50
      ELSEIF (LPIECE .GT. LWRAP+1) THEN
C
C       LPIECE SHOULD BE SET DOWN TO LWRAP.
C
         IDELTA = 0
         LPIECE = LWRAP
         DO 56 I=LPIECE+1,2,-1
            IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
               LPIECE = I-1
               IDELTA = 1
               GOTO 58
            ENDIF
   56    CONTINUE
   58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSE
C
C       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.
C       WE SHOULD DECREMENT LPIECE BY ONE.
C
         LPIECE = LPIECE - 1
         CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC  = NEXTC + LPIECE + 2
      ENDIF
C
C       PRINT
C
      DO 60 I=1,NUNIT
         WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
   60 CONTINUE
C
      IF (NEXTC .LE. LENMSG) GO TO 50
      RETURN
      END
*DECK XERSVE
      SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL,
     +   ICOUNT)
C***BEGIN PROLOGUE  XERSVE
C***SUBSIDIARY
C***PURPOSE  Record that an error has occurred.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3
C***TYPE      ALL (XERSVE-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C *Usage:
C
C        INTEGER  KFLAG, NERR, LEVEL, ICOUNT
C        CHARACTER * (len) LIBRAR, SUBROU, MESSG
C
C        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
C
C *Arguments:
C
C        LIBRAR :IN    is the library that the message is from.
C        SUBROU :IN    is the subroutine that the message is from.
C        MESSG  :IN    is the message to be saved.
C        KFLAG  :IN    indicates the action to be performed.
C                      when KFLAG > 0, the message in MESSG is saved.
C                      when KFLAG=0 the tables will be dumped and
C                      cleared.
C                      when KFLAG < 0, the tables will be dumped and
C                      not cleared.
C        NERR   :IN    is the error number.
C        LEVEL  :IN    is the error severity.
C        ICOUNT :OUT   the number of times this message has been seen,
C                      or zero if the table has overflowed and does not
C                      contain this message specifically.  When KFLAG=0,
C                      ICOUNT will not be altered.
C
C *Description:
C
C   Record that this error occurred and possibly dump and clear the
C   tables.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  I1MACH, XGETUA
C***REVISION HISTORY  (YYMMDD)
C   800319  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900413  Routine modified to remove reference to KFLAG.  (WRB)
C   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
C           sequence, use IF-THEN-ELSE, make number of saved entries
C           easily changeable, changed routine name from XERSAV to
C           XERSVE.  (RWC)
C   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERSVE
      PARAMETER (LENTAB=10)
      INTEGER LUN(5)
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
      CHARACTER*20 MESTAB(LENTAB), MES
      DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)
      SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
      DATA KOUNTX/0/, NMSG/0/
C***FIRST EXECUTABLE STATEMENT  XERSVE
C
      IF (KFLAG.LE.0) THEN
C
C        Dump the table.
C
         IF (NMSG.EQ.0) RETURN
C
C        Print to each unit.
C
         CALL XGETUA (LUN, NUNIT)
         DO 20 KUNIT = 1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C
C           Print the table header.
C
            WRITE (IUNIT,9000)
C
C           Print body of table.
C
            DO 10 I = 1,NMSG
               WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I),
     *            NERTAB(I),LEVTAB(I),KOUNT(I)
   10       CONTINUE
C
C           Print number of other errors.
C
            IF (KOUNTX.NE.0) WRITE (IUNIT,9020) KOUNTX
            WRITE (IUNIT,9030)
   20    CONTINUE
C
C        Clear the error tables.
C
         IF (KFLAG.EQ.0) THEN
            NMSG = 0
            KOUNTX = 0
         ENDIF
      ELSE
C
C        PROCESS A MESSAGE...
C        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
C
         LIB = LIBRAR
         SUB = SUBROU
         MES = MESSG
         DO 30 I = 1,NMSG
            IF (LIB.EQ.LIBTAB(I) .AND. SUB.EQ.SUBTAB(I) .AND.
     *         MES.EQ.MESTAB(I) .AND. NERR.EQ.NERTAB(I) .AND.
     *         LEVEL.EQ.LEVTAB(I)) THEN
                  KOUNT(I) = KOUNT(I) + 1
                  ICOUNT = KOUNT(I)
                  RETURN
            ENDIF
   30    CONTINUE
C
         IF (NMSG.LT.LENTAB) THEN
C
C           Empty slot found for new message.
C
            NMSG = NMSG + 1
            LIBTAB(I) = LIB
            SUBTAB(I) = SUB
            MESTAB(I) = MES
            NERTAB(I) = NERR
            LEVTAB(I) = LEVEL
            KOUNT (I) = 1
            ICOUNT    = 1
         ELSE
C
C           Table is full.
C
            KOUNTX = KOUNTX+1
            ICOUNT = 0
         ENDIF
      ENDIF
      RETURN
C
C     Formats.
C
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' /
     +   ' LIBRARY    SUBROUTINE MESSAGE START             NERR',
     +   '     LEVEL     COUNT')
 9010 FORMAT (1X,A,3X,A,3X,A,3I10)
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
 9030 FORMAT (1X)
      END
*DECK XGETUA
      SUBROUTINE XGETUA (IUNITA, N)
C***BEGIN PROLOGUE  XGETUA
C***PURPOSE  Return unit number(s) to which error messages are being
C            sent.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XGETUA-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        XGETUA may be called to determine the unit number or numbers
C        to which error messages are being sent.
C        These unit numbers may have been set by a call to XSETUN,
C        or a call to XSETUA, or may be a default value.
C
C     Description of Parameters
C      --Output--
C        IUNIT - an array of one to five unit numbers, depending
C                on the value of N.  A value of zero refers to the
C                default unit, as defined by the I1MACH machine
C                constant routine.  Only IUNIT(1),...,IUNIT(N) are
C                defined by XGETUA.  The values of IUNIT(N+1),...,
C                IUNIT(5) are not defined (for N .LT. 5) or altered
C                in any way by XGETUA.
C        N     - the number of units to which copies of the
C                error messages are being sent.  N will be in the
C                range from 1 to 5.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  J4SAVE
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XGETUA
      DIMENSION IUNITA(5)
C***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
*DECK ZABS
      DOUBLE PRECISION FUNCTION ZABS (ZR, ZI)
C***BEGIN PROLOGUE  ZABS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
C            ZBIRY
C***LIBRARY   SLATEC
C***TYPE      ALL (ZABS-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZABS COMPUTES THE ABSOLUTE VALUE OR MAGNITUDE OF A DOUBLE
C     PRECISION COMPLEX VARIABLE CMPLX(ZR,ZI)
C
C***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZABS
      DOUBLE PRECISION ZR, ZI, U, V, Q, S
C***FIRST EXECUTABLE STATEMENT  ZABS
      U = ABS(ZR)
      V = ABS(ZI)
      S = U + V
C-----------------------------------------------------------------------
C     S*1.0D0 MAKES AN UNNORMALIZED UNDERFLOW ON CDC MACHINES INTO A
C     TRUE FLOATING ZERO
C-----------------------------------------------------------------------
      S = S*1.0D+0
      IF (S.EQ.0.0D+0) GO TO 20
      IF (U.GT.V) GO TO 10
      Q = U/V
      ZABS = V*SQRT(1.D+0+Q*Q)
      RETURN
   10 Q = V/U
      ZABS = U*SQRT(1.D+0+Q*Q)
      RETURN
   20 ZABS = 0.0D+0
      RETURN
      END
*DECK ZACAI
      SUBROUTINE ZACAI (ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, RL, TOL,
     +   ELIM, ALIM)
C***BEGIN PROLOGUE  ZACAI
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZAIRY
C***LIBRARY   SLATEC
C***TYPE      ALL (CACAI-A, ZACAI-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
C
C         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
C                 MP=PI*MR*CMPLX(0.0,1.0)
C
C     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
C     HALF Z PLANE FOR USE WITH ZAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
C     ZACAI IS THE SAME AS ZACON WITH THE PARTS FOR LARGER ORDERS AND
C     RECURRENCE REMOVED. A RECURSIVE CALL TO ZACON CAN RESULT IF ZACON
C     IS CALLED FROM ZAIRY.
C
C***SEE ALSO  ZAIRY
C***ROUTINES CALLED  D1MACH, ZABS, ZASYI, ZBKNU, ZMLRI, ZS1S2, ZSERI
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZACAI
C     COMPLEX CSGN,CSPN,C1,C2,Y,Z,ZN,CY
      DOUBLE PRECISION ALIM, ARG, ASCLE, AZ, CSGNR, CSGNI, CSPNR,
     * CSPNI, C1R, C1I, C2R, C2I, CYR, CYI, DFNU, ELIM, FMR, FNU, PI,
     * RL, SGN, TOL, YY, YR, YI, ZR, ZI, ZNR, ZNI, D1MACH, ZABS
      INTEGER INU, IUF, KODE, MR, N, NN, NW, NZ
      DIMENSION YR(N), YI(N), CYR(2), CYI(2)
      EXTERNAL ZABS
      DATA PI / 3.14159265358979324D0 /
C***FIRST EXECUTABLE STATEMENT  ZACAI
      NZ = 0
      ZNR = -ZR
      ZNI = -ZI
      AZ = ZABS(ZR,ZI)
      NN = N
      DFNU = FNU + (N-1)
      IF (AZ.LE.2.0D0) GO TO 10
      IF (AZ*AZ*0.25D0.GT.DFNU+1.0D0) GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     POWER SERIES FOR THE I FUNCTION
C-----------------------------------------------------------------------
      CALL ZSERI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL, ELIM, ALIM)
      GO TO 40
   20 CONTINUE
      IF (AZ.LT.RL) GO TO 30
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
C-----------------------------------------------------------------------
      CALL ZASYI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, RL, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 80
      GO TO 40
   30 CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
C-----------------------------------------------------------------------
      CALL ZMLRI(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, TOL)
      IF(NW.LT.0) GO TO 80
   40 CONTINUE
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
C-----------------------------------------------------------------------
      CALL ZBKNU(ZNR, ZNI, FNU, KODE, 1, CYR, CYI, NW, TOL, ELIM, ALIM)
      IF (NW.NE.0) GO TO 80
      FMR = MR
      SGN = -DSIGN(PI,FMR)
      CSGNR = 0.0D0
      CSGNI = SGN
      IF (KODE.EQ.1) GO TO 50
      YY = -ZNI
      CSGNR = -CSGNI*SIN(YY)
      CSGNI = CSGNI*COS(YY)
   50 CONTINUE
C-----------------------------------------------------------------------
C     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = FNU
      ARG = (FNU-INU)*SGN
      CSPNR = COS(ARG)
      CSPNI = SIN(ARG)
      IF (MOD(INU,2).EQ.0) GO TO 60
      CSPNR = -CSPNR
      CSPNI = -CSPNI
   60 CONTINUE
      C1R = CYR(1)
      C1I = CYI(1)
      C2R = YR(1)
      C2I = YI(1)
      IF (KODE.EQ.1) GO TO 70
      IUF = 0
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
   70 CONTINUE
      YR(1) = CSPNR*C1R - CSPNI*C1I + CSGNR*C2R - CSGNI*C2I
      YI(1) = CSPNR*C1I + CSPNI*C1R + CSGNR*C2I + CSGNI*C2R
      RETURN
   80 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END
*DECK ZACON
      SUBROUTINE ZACON (ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, RL, FNUL,
     +   TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  ZACON
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CACON-A, ZACON-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZACON APPLIES THE ANALYTIC CONTINUATION FORMULA
C
C         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
C                 MP=PI*MR*CMPLX(0.0,1.0)
C
C     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
C     HALF Z PLANE
C
C***SEE ALSO  ZBESH, ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZBINU, ZBKNU, ZMLT, ZS1S2
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZACON
C     COMPLEX CK,CONE,CSCL,CSCR,CSGN,CSPN,CY,CZERO,C1,C2,RZ,SC1,SC2,ST,
C    *S1,S2,Y,Z,ZN
      DOUBLE PRECISION ALIM, ARG, ASCLE, AS2, AZN, BRY, BSCLE, CKI,
     * CKR, CONER, CPN, CSCL, CSCR, CSGNI, CSGNR, CSPNI, CSPNR,
     * CSR, CSRR, CSSR, CYI, CYR, C1I, C1M, C1R, C2I, C2R, ELIM, FMR,
     * FN, FNU, FNUL, PI, PTI, PTR, RAZN, RL, RZI, RZR, SC1I, SC1R,
     * SC2I, SC2R, SGN, SPN, STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR,
     * YY, ZEROR, ZI, ZNI, ZNR, ZR, D1MACH, ZABS
      INTEGER I, INU, IUF, KFLAG, KODE, MR, N, NN, NW, NZ
      DIMENSION YR(N), YI(N), CYR(2), CYI(2), CSSR(3), CSRR(3), BRY(3)
      EXTERNAL ZABS
      DATA PI / 3.14159265358979324D0 /
      DATA ZEROR,CONER / 0.0D0,1.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZACON
      NZ = 0
      ZNR = -ZR
      ZNI = -ZI
      NN = N
      CALL ZBINU(ZNR, ZNI, FNU, KODE, NN, YR, YI, NW, RL, FNUL, TOL,
     * ELIM, ALIM)
      IF (NW.LT.0) GO TO 90
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
C-----------------------------------------------------------------------
      NN = MIN(2,N)
      CALL ZBKNU(ZNR, ZNI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM)
      IF (NW.NE.0) GO TO 90
      S1R = CYR(1)
      S1I = CYI(1)
      FMR = MR
      SGN = -DSIGN(PI,FMR)
      CSGNR = ZEROR
      CSGNI = SGN
      IF (KODE.EQ.1) GO TO 10
      YY = -ZNI
      CPN = COS(YY)
      SPN = SIN(YY)
      CALL ZMLT(CSGNR, CSGNI, CPN, SPN, CSGNR, CSGNI)
   10 CONTINUE
C-----------------------------------------------------------------------
C     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = FNU
      ARG = (FNU-INU)*SGN
      CPN = COS(ARG)
      SPN = SIN(ARG)
      CSPNR = CPN
      CSPNI = SPN
      IF (MOD(INU,2).EQ.0) GO TO 20
      CSPNR = -CSPNR
      CSPNI = -CSPNI
   20 CONTINUE
      IUF = 0
      C1R = S1R
      C1I = S1I
      C2R = YR(1)
      C2I = YI(1)
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      IF (KODE.EQ.1) GO TO 30
      CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
      SC1R = C1R
      SC1I = C1I
   30 CONTINUE
      CALL ZMLT(CSPNR, CSPNI, C1R, C1I, STR, STI)
      CALL ZMLT(CSGNR, CSGNI, C2R, C2I, PTR, PTI)
      YR(1) = STR + PTR
      YI(1) = STI + PTI
      IF (N.EQ.1) RETURN
      CSPNR = -CSPNR
      CSPNI = -CSPNI
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = S2R
      C1I = S2I
      C2R = YR(2)
      C2I = YI(2)
      IF (KODE.EQ.1) GO TO 40
      CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
      NZ = NZ + NW
      SC2R = C1R
      SC2I = C1I
   40 CONTINUE
      CALL ZMLT(CSPNR, CSPNI, C1R, C1I, STR, STI)
      CALL ZMLT(CSGNR, CSGNI, C2R, C2I, PTR, PTI)
      YR(2) = STR + PTR
      YI(2) = STI + PTI
      IF (N.EQ.2) RETURN
      CSPNR = -CSPNR
      CSPNI = -CSPNI
      AZN = ZABS(ZNR,ZNI)
      RAZN = 1.0D0/AZN
      STR = ZNR*RAZN
      STI = -ZNI*RAZN
      RZR = (STR+STR)*RAZN
      RZI = (STI+STI)*RAZN
      FN = FNU + 1.0D0
      CKR = FN*RZR
      CKI = FN*RZI
C-----------------------------------------------------------------------
C     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
C-----------------------------------------------------------------------
      CSCL = 1.0D0/TOL
      CSCR = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CSCR
      CSRR(1) = CSCR
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = ASCLE
      BRY(2) = 1.0D0/ASCLE
      BRY(3) = D1MACH(2)
      AS2 = ZABS(S2R,S2I)
      KFLAG = 2
      IF (AS2.GT.BRY(1)) GO TO 50
      KFLAG = 1
      GO TO 60
   50 CONTINUE
      IF (AS2.LT.BRY(2)) GO TO 60
      KFLAG = 3
   60 CONTINUE
      BSCLE = BRY(KFLAG)
      S1R = S1R*CSSR(KFLAG)
      S1I = S1I*CSSR(KFLAG)
      S2R = S2R*CSSR(KFLAG)
      S2I = S2I*CSSR(KFLAG)
      CSR = CSRR(KFLAG)
      DO 80 I=3,N
        STR = S2R
        STI = S2I
        S2R = CKR*STR - CKI*STI + S1R
        S2I = CKR*STI + CKI*STR + S1I
        S1R = STR
        S1I = STI
        C1R = S2R*CSR
        C1I = S2I*CSR
        STR = C1R
        STI = C1I
        C2R = YR(I)
        C2I = YI(I)
        IF (KODE.EQ.1) GO TO 70
        IF (IUF.LT.0) GO TO 70
        CALL ZS1S2(ZNR, ZNI, C1R, C1I, C2R, C2I, NW, ASCLE, ALIM, IUF)
        NZ = NZ + NW
        SC1R = SC2R
        SC1I = SC2I
        SC2R = C1R
        SC2I = C1I
        IF (IUF.NE.3) GO TO 70
        IUF = -4
        S1R = SC1R*CSSR(KFLAG)
        S1I = SC1I*CSSR(KFLAG)
        S2R = SC2R*CSSR(KFLAG)
        S2I = SC2I*CSSR(KFLAG)
        STR = SC2R
        STI = SC2I
   70   CONTINUE
        PTR = CSPNR*C1R - CSPNI*C1I
        PTI = CSPNR*C1I + CSPNI*C1R
        YR(I) = PTR + CSGNR*C2R - CSGNI*C2I
        YI(I) = PTI + CSGNR*C2I + CSGNI*C2R
        CKR = CKR + RZR
        CKI = CKI + RZI
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        IF (KFLAG.GE.3) GO TO 80
        PTR = ABS(C1R)
        PTI = ABS(C1I)
        C1M = MAX(PTR,PTI)
        IF (C1M.LE.BSCLE) GO TO 80
        KFLAG = KFLAG + 1
        BSCLE = BRY(KFLAG)
        S1R = S1R*CSR
        S1I = S1I*CSR
        S2R = STR
        S2I = STI
        S1R = S1R*CSSR(KFLAG)
        S1I = S1I*CSSR(KFLAG)
        S2R = S2R*CSSR(KFLAG)
        S2I = S2I*CSSR(KFLAG)
        CSR = CSRR(KFLAG)
   80 CONTINUE
      RETURN
   90 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END
*DECK ZAIRY
      SUBROUTINE ZAIRY (ZR, ZI, ID, KODE, AIR, AII, NZ, IERR)
C***BEGIN PROLOGUE  ZAIRY
C***PURPOSE  Compute the Airy function Ai(z) or its derivative dAi/dz
C            for complex argument z.  A scaling option is available
C            to help avoid underflow and overflow.
C***LIBRARY   SLATEC
C***CATEGORY  C10D
C***TYPE      COMPLEX (CAIRY-C, ZAIRY-C)
C***KEYWORDS  AIRY FUNCTION, BESSEL FUNCTION OF ORDER ONE THIRD,
C             BESSEL FUNCTION OF ORDER TWO THIRDS
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C         On KODE=1, ZAIRY computes the complex Airy function Ai(z)
C         or its derivative dAi/dz on ID=0 or ID=1 respectively. On
C         KODE=2, a scaling option exp(zeta)*Ai(z) or exp(zeta)*dAi/dz
C         is provided to remove the exponential decay in -pi/3<arg(z)
C         <pi/3 and the exponential growth in pi/3<abs(arg(z))<pi where
C         zeta=(2/3)*z**(3/2).
C
C         While the Airy functions Ai(z) and dAi/dz are analytic in
C         the whole z-plane, the corresponding scaled functions defined
C         for KODE=2 have a cut along the negative real axis.
C
C         Input
C           ZR     - DOUBLE PRECISION real part of argument Z
C           ZI     - DOUBLE PRECISION imag part of argument Z
C           ID     - Order of derivative, ID=0 or ID=1
C           KODE   - A parameter to indicate the scaling option
C                    KODE=1  returns
C                            AI=Ai(z)  on ID=0
C                            AI=dAi/dz on ID=1
C                            at z=Z
C                        =2  returns
C                            AI=exp(zeta)*Ai(z)  on ID=0
C                            AI=exp(zeta)*dAi/dz on ID=1
C                            at z=Z where zeta=(2/3)*z**(3/2)
C
C         Output
C           AIR    - DOUBLE PRECISION real part of result
C           AII    - DOUBLE PRECISION imag part of result
C           NZ     - Underflow indicator
C                    NZ=0    Normal return
C                    NZ=1    AI=0 due to underflow in
C                            -pi/3<arg(Z)<pi/3 on KODE=1
C           IERR   - Error flag
C                    IERR=0  Normal return     - COMPUTATION COMPLETED
C                    IERR=1  Input error       - NO COMPUTATION
C                    IERR=2  Overflow          - NO COMPUTATION
C                            (Re(Z) too large with KODE=1)
C                    IERR=3  Precision warning - COMPUTATION COMPLETED
C                            (Result has less than half precision)
C                    IERR=4  Precision error   - NO COMPUTATION
C                            (Result has no precision)
C                    IERR=5  Algorithmic error - NO COMPUTATION
C                            (Termination condition not met)
C
C *Long Description:
C
C         Ai(z) and dAi/dz are computed from K Bessel functions by
C
C                Ai(z) =  c*sqrt(z)*K(1/3,zeta)
C               dAi/dz = -c*   z   *K(2/3,zeta)
C                    c =  1/(pi*sqrt(3))
C                 zeta =  (2/3)*z**(3/2)
C
C         when abs(z)>1 and from power series when abs(z)<=1.
C
C         In most complex variable computation, one must evaluate ele-
C         mentary functions.  When the magnitude of Z is large, losses
C         of significance by argument reduction occur.  Consequently, if
C         the magnitude of ZETA=(2/3)*Z**(3/2) exceeds U1=SQRT(0.5/UR),
C         then losses exceeding half precision are likely and an error
C         flag IERR=3 is triggered where UR=MAX(D1MACH(4),1.0D-18) is
C         double precision unit roundoff limited to 18 digits precision.
C         Also, if the magnitude of ZETA is larger than U2=0.5/UR, then
C         all significance is lost and IERR=4.  In order to use the INT
C         function, ZETA must be further restricted not to exceed
C         U3=I1MACH(9)=LARGEST INTEGER.  Thus, the magnitude of ZETA
C         must be restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2,
C         and U3 are approximately 2.0E+3, 4.2E+6, 2.1E+9 in single
C         precision and 4.7E+7, 2.3E+15, 2.1E+9 in double precision.
C         This makes U2 limiting is single precision and U3 limiting
C         in double precision.  This means that the magnitude of Z
C         cannot exceed approximately 3.4E+4 in single precision and
C         2.1E+6 in double precision.  This also means that one can
C         expect to retain, in the worst cases on 32-bit machines,
C         no digits in single precision and only 6 digits in double
C         precision.
C
C         The approximate relative error in the magnitude of a complex
C         Bessel function can be expressed as P*10**S where P=MAX(UNIT
C         ROUNDOFF,1.0E-18) is the nominal precision and 10**S repre-
C         sents the increase in error due to argument reduction in the
C         elementary functions.  Here, S=MAX(1,ABS(LOG10(ABS(Z))),
C         ABS(LOG10(FNU))) approximately (i.e., S=MAX(1,ABS(EXPONENT OF
C         ABS(Z),ABS(EXPONENT OF FNU)) ).  However, the phase angle may
C         have only absolute accuracy.  This is most likely to occur
C         when one component (in magnitude) is larger than the other by
C         several orders of magnitude.  If one component is 10**K larger
C         than the other, then one can expect only MAX(ABS(LOG10(P))-K,
C         0) significant digits; or, stated another way, when K exceeds
C         the exponent of P, no significant digits remain in the smaller
C         component.  However, the phase angle retains absolute accuracy
C         because, in complex arithmetic with precision P, the smaller
C         component will not (as a rule) decrease below P times the
C         magnitude of the larger component. In these extreme cases,
C         the principal phase angle is on the order of +P, -P, PI/2-P,
C         or -PI/2+P.
C
C***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe-
C                 matical Functions, National Bureau of Standards
C                 Applied Mathematics Series 55, U. S. Department
C                 of Commerce, Tenth Printing (1972) or later.
C               2. D. E. Amos, Computation of Bessel Functions of
C                 Complex Argument and Large Order, Report SAND83-0643,
C                 Sandia National Laboratories, Albuquerque, NM, May
C                 1983.
C               3. D. E. Amos, A Subroutine Package for Bessel Functions
C                 of a Complex Argument and Nonnegative Order, Report
C                 SAND85-1018, Sandia National Laboratory, Albuquerque,
C                 NM, May 1985.
C               4. D. E. Amos, A portable package for Bessel functions
C                 of a complex argument and nonnegative order, ACM
C                 Transactions on Mathematical Software, 12 (September
C                 1986), pp. 265-273.
C
C***ROUTINES CALLED  D1MACH, I1MACH, ZABS, ZACAI, ZBKNU, ZEXP, ZSQRT
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   890801  REVISION DATE from Version 3.2
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   920128  Category corrected.  (WRB)
C   920811  Prologue revised.  (DWL)
C   930122  Added ZEXP and ZSQRT to EXTERNAL statement.  (RWC)
C***END PROLOGUE  ZAIRY
C     COMPLEX AI,CONE,CSQ,CY,S1,S2,TRM1,TRM2,Z,ZTA,Z3
      DOUBLE PRECISION AA, AD, AII, AIR, AK, ALIM, ATRM, AZ, AZ3, BK,
     * CC, CK, COEF, CONEI, CONER, CSQI, CSQR, CYI, CYR, C1, C2, DIG,
     * DK, D1, D2, ELIM, FID, FNU, PTR, RL, R1M5, SFAC, STI, STR,
     * S1I, S1R, S2I, S2R, TOL, TRM1I, TRM1R, TRM2I, TRM2R, TTH, ZEROI,
     * ZEROR, ZI, ZR, ZTAI, ZTAR, Z3I, Z3R, D1MACH, ZABS, ALAZ, BB
      INTEGER ID, IERR, IFLAG, K, KODE, K1, K2, MR, NN, NZ, I1MACH
      DIMENSION CYR(1), CYI(1)
      EXTERNAL ZABS, ZEXP, ZSQRT
      DATA TTH, C1, C2, COEF /6.66666666666666667D-01,
     * 3.55028053887817240D-01,2.58819403792806799D-01,
     * 1.83776298473930683D-01/
      DATA ZEROR, ZEROI, CONER, CONEI /0.0D0,0.0D0,1.0D0,0.0D0/
C***FIRST EXECUTABLE STATEMENT  ZAIRY
      IERR = 0
      NZ=0
      IF (ID.LT.0 .OR. ID.GT.1) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (IERR.NE.0) RETURN
      AZ = ZABS(ZR,ZI)
      TOL = MAX(D1MACH(4),1.0D-18)
      FID = ID
      IF (AZ.GT.1.0D0) GO TO 70
C-----------------------------------------------------------------------
C     POWER SERIES FOR ABS(Z).LE.1.
C-----------------------------------------------------------------------
      S1R = CONER
      S1I = CONEI
      S2R = CONER
      S2I = CONEI
      IF (AZ.LT.TOL) GO TO 170
      AA = AZ*AZ
      IF (AA.LT.TOL/AZ) GO TO 40
      TRM1R = CONER
      TRM1I = CONEI
      TRM2R = CONER
      TRM2I = CONEI
      ATRM = 1.0D0
      STR = ZR*ZR - ZI*ZI
      STI = ZR*ZI + ZI*ZR
      Z3R = STR*ZR - STI*ZI
      Z3I = STR*ZI + STI*ZR
      AZ3 = AZ*AA
      AK = 2.0D0 + FID
      BK = 3.0D0 - FID - FID
      CK = 4.0D0 - FID
      DK = 3.0D0 + FID + FID
      D1 = AK*DK
      D2 = BK*CK
      AD = MIN(D1,D2)
      AK = 24.0D0 + 9.0D0*FID
      BK = 30.0D0 - 9.0D0*FID
      DO 30 K=1,25
        STR = (TRM1R*Z3R-TRM1I*Z3I)/D1
        TRM1I = (TRM1R*Z3I+TRM1I*Z3R)/D1
        TRM1R = STR
        S1R = S1R + TRM1R
        S1I = S1I + TRM1I
        STR = (TRM2R*Z3R-TRM2I*Z3I)/D2
        TRM2I = (TRM2R*Z3I+TRM2I*Z3R)/D2
        TRM2R = STR
        S2R = S2R + TRM2R
        S2I = S2I + TRM2I
        ATRM = ATRM*AZ3/AD
        D1 = D1 + AK
        D2 = D2 + BK
        AD = MIN(D1,D2)
        IF (ATRM.LT.TOL*AD) GO TO 40
        AK = AK + 18.0D0
        BK = BK + 18.0D0
   30 CONTINUE
   40 CONTINUE
      IF (ID.EQ.1) GO TO 50
      AIR = S1R*C1 - C2*(ZR*S2R-ZI*S2I)
      AII = S1I*C1 - C2*(ZR*S2I+ZI*S2R)
      IF (KODE.EQ.1) RETURN
      CALL ZSQRT(ZR, ZI, STR, STI)
      ZTAR = TTH*(ZR*STR-ZI*STI)
      ZTAI = TTH*(ZR*STI+ZI*STR)
      CALL ZEXP(ZTAR, ZTAI, STR, STI)
      PTR = AIR*STR - AII*STI
      AII = AIR*STI + AII*STR
      AIR = PTR
      RETURN
   50 CONTINUE
      AIR = -S2R*C2
      AII = -S2I*C2
      IF (AZ.LE.TOL) GO TO 60
      STR = ZR*S1R - ZI*S1I
      STI = ZR*S1I + ZI*S1R
      CC = C1/(1.0D0+FID)
      AIR = AIR + CC*(STR*ZR-STI*ZI)
      AII = AII + CC*(STR*ZI+STI*ZR)
   60 CONTINUE
      IF (KODE.EQ.1) RETURN
      CALL ZSQRT(ZR, ZI, STR, STI)
      ZTAR = TTH*(ZR*STR-ZI*STI)
      ZTAI = TTH*(ZR*STI+ZI*STR)
      CALL ZEXP(ZTAR, ZTAI, STR, STI)
      PTR = STR*AIR - STI*AII
      AII = STR*AII + STI*AIR
      AIR = PTR
      RETURN
C-----------------------------------------------------------------------
C     CASE FOR ABS(Z).GT.1.0
C-----------------------------------------------------------------------
   70 CONTINUE
      FNU = (1.0D0+FID)/3.0D0
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0D-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C-----------------------------------------------------------------------
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN(ABS(K1),ABS(K2))
      ELIM = 2.303D0*(K*R1M5-3.0D0)
      K1 = I1MACH(14) - 1
      AA = R1M5*K1
      DIG = MIN(AA,18.0D0)
      AA = AA*2.303D0
      ALIM = ELIM + MAX(-AA,-41.45D0)
      RL = 1.2D0*DIG + 3.0D0
      ALAZ = LOG(AZ)
C-----------------------------------------------------------------------
C     TEST FOR PROPER RANGE
C-----------------------------------------------------------------------
      AA=0.5D0/TOL
      BB=I1MACH(9)*0.5D0
      AA=MIN(AA,BB)
      AA=AA**TTH
      IF (AZ.GT.AA) GO TO 260
      AA=SQRT(AA)
      IF (AZ.GT.AA) IERR=3
      CALL ZSQRT(ZR, ZI, CSQR, CSQI)
      ZTAR = TTH*(ZR*CSQR-ZI*CSQI)
      ZTAI = TTH*(ZR*CSQI+ZI*CSQR)
C-----------------------------------------------------------------------
C     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL
C-----------------------------------------------------------------------
      IFLAG = 0
      SFAC = 1.0D0
      AK = ZTAI
      IF (ZR.GE.0.0D0) GO TO 80
      BK = ZTAR
      CK = -ABS(BK)
      ZTAR = CK
      ZTAI = AK
   80 CONTINUE
      IF (ZI.NE.0.0D0) GO TO 90
      IF (ZR.GT.0.0D0) GO TO 90
      ZTAR = 0.0D0
      ZTAI = AK
   90 CONTINUE
      AA = ZTAR
      IF (AA.GE.0.0D0 .AND. ZR.GT.0.0D0) GO TO 110
      IF (KODE.EQ.2) GO TO 100
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
      IF (AA.GT.(-ALIM)) GO TO 100
      AA = -AA + 0.25D0*ALAZ
      IFLAG = 1
      SFAC = TOL
      IF (AA.GT.ELIM) GO TO 270
  100 CONTINUE
C-----------------------------------------------------------------------
C     CBKNU AND CACON RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
C-----------------------------------------------------------------------
      MR = 1
      IF (ZI.LT.0.0D0) MR = -1
      CALL ZACAI(ZTAR, ZTAI, FNU, KODE, MR, 1, CYR, CYI, NN, RL, TOL,
     * ELIM, ALIM)
      IF (NN.LT.0) GO TO 280
      NZ = NZ + NN
      GO TO 130
  110 CONTINUE
      IF (KODE.EQ.2) GO TO 120
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
      IF (AA.LT.ALIM) GO TO 120
      AA = -AA - 0.25D0*ALAZ
      IFLAG = 2
      SFAC = 1.0D0/TOL
      IF (AA.LT.(-ELIM)) GO TO 210
  120 CONTINUE
      CALL ZBKNU(ZTAR, ZTAI, FNU, KODE, 1, CYR, CYI, NZ, TOL, ELIM,
     * ALIM)
  130 CONTINUE
      S1R = CYR(1)*COEF
      S1I = CYI(1)*COEF
      IF (IFLAG.NE.0) GO TO 150
      IF (ID.EQ.1) GO TO 140
      AIR = CSQR*S1R - CSQI*S1I
      AII = CSQR*S1I + CSQI*S1R
      RETURN
  140 CONTINUE
      AIR = -(ZR*S1R-ZI*S1I)
      AII = -(ZR*S1I+ZI*S1R)
      RETURN
  150 CONTINUE
      S1R = S1R*SFAC
      S1I = S1I*SFAC
      IF (ID.EQ.1) GO TO 160
      STR = S1R*CSQR - S1I*CSQI
      S1I = S1R*CSQI + S1I*CSQR
      S1R = STR
      AIR = S1R/SFAC
      AII = S1I/SFAC
      RETURN
  160 CONTINUE
      STR = -(S1R*ZR-S1I*ZI)
      S1I = -(S1R*ZI+S1I*ZR)
      S1R = STR
      AIR = S1R/SFAC
      AII = S1I/SFAC
      RETURN
  170 CONTINUE
      AA = 1.0D+3*D1MACH(1)
      S1R = ZEROR
      S1I = ZEROI
      IF (ID.EQ.1) GO TO 190
      IF (AZ.LE.AA) GO TO 180
      S1R = C2*ZR
      S1I = C2*ZI
  180 CONTINUE
      AIR = C1 - S1R
      AII = -S1I
      RETURN
  190 CONTINUE
      AIR = -C2
      AII = 0.0D0
      AA = SQRT(AA)
      IF (AZ.LE.AA) GO TO 200
      S1R = 0.5D0*(ZR*ZR-ZI*ZI)
      S1I = ZR*ZI
  200 CONTINUE
      AIR = AIR + C1*S1R
      AII = AII + C1*S1I
      RETURN
  210 CONTINUE
      NZ = 1
      AIR = ZEROR
      AII = ZEROI
      RETURN
  270 CONTINUE
      NZ = 0
      IERR=2
      RETURN
  280 CONTINUE
      IF(NN.EQ.(-1)) GO TO 270
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      IERR=4
      NZ=0
      RETURN
      END
*DECK ZASYI
      SUBROUTINE ZASYI (ZR, ZI, FNU, KODE, N, YR, YI, NZ, RL, TOL, ELIM,
     +   ALIM)
C***BEGIN PROLOGUE  ZASYI
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CASYI-A, ZASYI-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
C     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE ABS(Z) IN THE
C     REGION ABS(Z).GT.MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN.
C     NZ.LT.0 INDICATES AN OVERFLOW ON KODE=1.
C
C***SEE ALSO  ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZDIV, ZEXP, ZMLT, ZSQRT
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   930122  Added ZEXP and ZSQRT to EXTERNAL statement.  (RWC)
C***END PROLOGUE  ZASYI
C     COMPLEX AK1,CK,CONE,CS1,CS2,CZ,CZERO,DK,EZ,P1,RZ,S2,Y,Z
      DOUBLE PRECISION AA, AEZ, AK, AK1I, AK1R, ALIM, ARG, ARM, ATOL,
     * AZ, BB, BK, CKI, CKR, CONEI, CONER, CS1I, CS1R, CS2I, CS2R, CZI,
     * CZR, DFNU, DKI, DKR, DNU2, ELIM, EZI, EZR, FDN, FNU, PI, P1I,
     * P1R, RAZ, RL, RTPI, RTR1, RZI, RZR, S, SGN, SQK, STI, STR, S2I,
     * S2R, TOL, TZI, TZR, YI, YR, ZEROI, ZEROR, ZI, ZR, D1MACH, ZABS
      INTEGER I, IB, IL, INU, J, JL, K, KODE, KODED, M, N, NN, NZ
      DIMENSION YR(N), YI(N)
      EXTERNAL ZABS, ZEXP, ZSQRT
      DATA PI, RTPI  /3.14159265358979324D0 , 0.159154943091895336D0 /
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZASYI
      NZ = 0
      AZ = ZABS(ZR,ZI)
      ARM = 1.0D+3*D1MACH(1)
      RTR1 = SQRT(ARM)
      IL = MIN(2,N)
      DFNU = FNU + (N-IL)
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
      RAZ = 1.0D0/AZ
      STR = ZR*RAZ
      STI = -ZI*RAZ
      AK1R = RTPI*STR*RAZ
      AK1I = RTPI*STI*RAZ
      CALL ZSQRT(AK1R, AK1I, AK1R, AK1I)
      CZR = ZR
      CZI = ZI
      IF (KODE.NE.2) GO TO 10
      CZR = ZEROR
      CZI = ZI
   10 CONTINUE
      IF (ABS(CZR).GT.ELIM) GO TO 100
      DNU2 = DFNU + DFNU
      KODED = 1
      IF ((ABS(CZR).GT.ALIM) .AND. (N.GT.2)) GO TO 20
      KODED = 0
      CALL ZEXP(CZR, CZI, STR, STI)
      CALL ZMLT(AK1R, AK1I, STR, STI, AK1R, AK1I)
   20 CONTINUE
      FDN = 0.0D0
      IF (DNU2.GT.RTR1) FDN = DNU2*DNU2
      EZR = ZR*8.0D0
      EZI = ZI*8.0D0
C-----------------------------------------------------------------------
C     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE
C     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE
C     EXPANSION FOR THE IMAGINARY PART.
C-----------------------------------------------------------------------
      AEZ = 8.0D0*AZ
      S = TOL/AEZ
      JL = RL+RL + 2
      P1R = ZEROR
      P1I = ZEROI
      IF (ZI.EQ.0.0D0) GO TO 30
C-----------------------------------------------------------------------
C     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF
C     SIGNIFICANCE WHEN FNU OR N IS LARGE
C-----------------------------------------------------------------------
      INU = FNU
      ARG = (FNU-INU)*PI
      INU = INU + N - IL
      AK = -SIN(ARG)
      BK = COS(ARG)
      IF (ZI.LT.0.0D0) BK = -BK
      P1R = AK
      P1I = BK
      IF (MOD(INU,2).EQ.0) GO TO 30
      P1R = -P1R
      P1I = -P1I
   30 CONTINUE
      DO 70 K=1,IL
        SQK = FDN - 1.0D0
        ATOL = S*ABS(SQK)
        SGN = 1.0D0
        CS1R = CONER
        CS1I = CONEI
        CS2R = CONER
        CS2I = CONEI
        CKR = CONER
        CKI = CONEI
        AK = 0.0D0
        AA = 1.0D0
        BB = AEZ
        DKR = EZR
        DKI = EZI
        DO 40 J=1,JL
          CALL ZDIV(CKR, CKI, DKR, DKI, STR, STI)
          CKR = STR*SQK
          CKI = STI*SQK
          CS2R = CS2R + CKR
          CS2I = CS2I + CKI
          SGN = -SGN
          CS1R = CS1R + CKR*SGN
          CS1I = CS1I + CKI*SGN
          DKR = DKR + EZR
          DKI = DKI + EZI
          AA = AA*ABS(SQK)/BB
          BB = BB + AEZ
          AK = AK + 8.0D0
          SQK = SQK - AK
          IF (AA.LE.ATOL) GO TO 50
   40   CONTINUE
        GO TO 110
   50   CONTINUE
        S2R = CS1R
        S2I = CS1I
        IF (ZR+ZR.GE.ELIM) GO TO 60
        TZR = ZR + ZR
        TZI = ZI + ZI
        CALL ZEXP(-TZR, -TZI, STR, STI)
        CALL ZMLT(STR, STI, P1R, P1I, STR, STI)
        CALL ZMLT(STR, STI, CS2R, CS2I, STR, STI)
        S2R = S2R + STR
        S2I = S2I + STI
   60   CONTINUE
        FDN = FDN + 8.0D0*DFNU + 4.0D0
        P1R = -P1R
        P1I = -P1I
        M = N - IL + K
        YR(M) = S2R*AK1R - S2I*AK1I
        YI(M) = S2R*AK1I + S2I*AK1R
   70 CONTINUE
      IF (N.LE.2) RETURN
      NN = N
      K = NN - 2
      AK = K
      STR = ZR*RAZ
      STI = -ZI*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      IB = 3
      DO 80 I=IB,NN
        YR(K) = (AK+FNU)*(RZR*YR(K+1)-RZI*YI(K+1)) + YR(K+2)
        YI(K) = (AK+FNU)*(RZR*YI(K+1)+RZI*YR(K+1)) + YI(K+2)
        AK = AK - 1.0D0
        K = K - 1
   80 CONTINUE
      IF (KODED.EQ.0) RETURN
      CALL ZEXP(CZR, CZI, CKR, CKI)
      DO 90 I=1,NN
        STR = YR(I)*CKR - YI(I)*CKI
        YI(I) = YR(I)*CKI + YI(I)*CKR
        YR(I) = STR
   90 CONTINUE
      RETURN
  100 CONTINUE
      NZ = -1
      RETURN
  110 CONTINUE
      NZ=-2
      RETURN
      END
      SUBROUTINE ZBESH(ZR, ZI, FNU, KODE, M, N, CYR, CYI, NZ, IERR)
C***BEGIN PROLOGUE  ZBESH
C***DATE WRITTEN   830501   (YYMMDD)
C***REVISION DATE  890801   (YYMMDD)
C***CATEGORY NO.  B5K
C***KEYWORDS  H-BESSEL FUNCTIONS,BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
C             BESSEL FUNCTIONS OF THIRD KIND,HANKEL FUNCTIONS
C***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C***PURPOSE  TO COMPUTE THE H-BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C         ON KODE=1, ZBESH COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         HANKEL (BESSEL) FUNCTIONS CY(J)=H(M,FNU+J-1,Z) FOR KINDS M=1
C         OR 2, REAL, NONNEGATIVE ORDERS FNU+J-1, J=1,...,N, AND COMPLEX
C         Z.NE.CMPLX(0.0,0.0) IN THE CUT PLANE -PI.LT.ARG(Z).LE.PI.
C         ON KODE=2, ZBESH RETURNS THE SCALED HANKEL FUNCTIONS
C
C         CY(I)=EXP(-MM*Z*I)*H(M,FNU+J-1,Z)       MM=3-2*M,   I**2=-1.
C
C         WHICH REMOVES THE EXPONENTIAL BEHAVIOR IN BOTH THE UPPER AND
C         LOWER HALF PLANES. DEFINITIONS AND NOTATION ARE FOUND IN THE
C         NBS HANDBOOK OF MATHEMATICAL FUNCTIONS (REF. 1).
C
C         INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
C           ZR,ZI  - Z=CMPLX(ZR,ZI), Z.NE.CMPLX(0.0D0,0.0D0),
C                    -PT.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL H FUNCTION, FNU.GE.0.0D0
C           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
C                    KODE= 1  RETURNS
C                             CY(J)=H(M,FNU+J-1,Z),   J=1,...,N
C                        = 2  RETURNS
C                             CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))
C                                  J=1,...,N  ,  I**2=-1
C           M      - KIND OF HANKEL FUNCTION, M=1 OR 2
C           N      - NUMBER OF MEMBERS IN THE SEQUENCE, N.GE.1
C
C         OUTPUT     CYR,CYI ARE DOUBLE PRECISION
C           CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
C                    CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
C                    CY(J)=H(M,FNU+J-1,Z)  OR
C                    CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))  J=1,...,N
C                    DEPENDING ON KODE, I**2=-1.
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , FIRST NZ COMPONENTS OF CY SET TO ZERO DUE
C                              TO UNDERFLOW, CY(J)=CMPLX(0.0D0,0.0D0)
C                              J=1,...,NZ WHEN Y.GT.0.0 AND M=1 OR
C                              Y.LT.0.0 AND M=2. FOR THE COMPLMENTARY
C                              HALF PLANES, NZ STATES ONLY THE NUMBER
C                              OF UNDERFLOWS.
C           IERR   - ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION, FNU TOO
C                            LARGE OR CABS(Z) TOO SMALL OR BOTH
C                    IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C***LONG DESCRIPTION
C
C         THE COMPUTATION IS CARRIED OUT BY THE RELATION
C
C         H(M,FNU,Z)=(1/MP)*EXP(-MP*FNU)*K(FNU,Z*EXP(-MP))
C             MP=MM*HPI*I,  MM=3-2*M,  HPI=PI/2,  I**2=-1
C
C         FOR M=1 OR 2 WHERE THE K BESSEL FUNCTION IS COMPUTED FOR THE
C         RIGHT HALF PLANE RE(Z).GE.0.0. THE K FUNCTION IS CONTINUED
C         TO THE LEFT HALF PLANE BY THE RELATION
C
C         K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
C         MP=MR*PI*I, MR=+1 OR -1, RE(Z).GT.0, I**2=-1
C
C         WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.
C
C         EXPONENTIAL DECAY OF H(M,FNU,Z) OCCURS IN THE UPPER HALF Z
C         PLANE FOR M=1 AND THE LOWER HALF Z PLANE FOR M=2.  EXPONENTIAL
C         GROWTH OCCURS IN THE COMPLEMENTARY HALF PLANES.  SCALING
C         BY EXP(-MM*Z*I) REMOVES THE EXPONENTIAL BEHAVIOR IN THE
C         WHOLE Z PLANE FOR Z TO INFINITY.
C
C         FOR NEGATIVE ORDERS,THE FORMULAE
C
C               H(1,-FNU,Z) = H(1,FNU,Z)*CEXP( PI*FNU*I)
C               H(2,-FNU,Z) = H(2,FNU,Z)*CEXP(-PI*FNU*I)
C                         I**2=-1
C
C         CAN BE USED.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=3 IS TRIGGERED WHERE UR=DMAX1(D1MACH(4),1.0D-18) IS
C         DOUBLE PRECISION UNIT ROUNDOFF LIMITED TO 18 DIGITS PRECISION.
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=I1MACH(9). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0D-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
C                 MATH. SOFTWARE, 1986
C
C***ROUTINES CALLED  ZACON,ZBKNU,ZBUNK,ZUOIK,ZABS,I1MACH,D1MACH
C***END PROLOGUE  ZBESH
C
C     COMPLEX CY,Z,ZN,ZT,CSGN
      DOUBLE PRECISION AA, ALIM, ALN, ARG, AZ, CYI, CYR, DIG, ELIM,
     * FMM, FN, FNU, FNUL, HPI, RHPI, RL, R1M5, SGN, STR, TOL, UFL, ZI,
     * ZNI, ZNR, ZR, ZTI, D1MACH, ZABS, BB, ASCLE, RTOL, ATOL, STI,
     * CSGNR, CSGNI
      INTEGER I, IERR, INU, INUH, IR, K, KODE, K1, K2, M,
     * MM, MR, N, NN, NUF, NW, NZ, I1MACH
      DIMENSION CYR(N), CYI(N)

c 110312 (BTD) add EXTERNAL statement:
      EXTERNAL ZABS
c------------------------------------
C
      DATA HPI /1.57079632679489662D0/
C
C***FIRST EXECUTABLE STATEMENT  ZBESH
      IERR = 0
      NZ=0
      IF (ZR.EQ.0.0D0 .AND. ZI.EQ.0.0D0) IERR=1
      IF (FNU.LT.0.0D0) IERR=1
      IF (M.LT.1 .OR. M.GT.2) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
      NN = N
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
C-----------------------------------------------------------------------
      TOL = DMAX1(D1MACH(4),1.0D-18)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN0(IABS(K1),IABS(K2))
      ELIM = 2.303D0*(DBLE(FLOAT(K))*R1M5-3.0D0)
      K1 = I1MACH(14) - 1
      AA = R1M5*DBLE(FLOAT(K1))
      DIG = DMIN1(AA,18.0D0)
      AA = AA*2.303D0
      ALIM = ELIM + DMAX1(-AA,-41.45D0)
      FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
      RL = 1.2D0*DIG + 3.0D0
      FN = FNU + DBLE(FLOAT(NN-1))
      MM = 3 - M - M
      FMM = DBLE(FLOAT(MM))
      ZNR = FMM*ZI
      ZNI = -FMM*ZR
C-----------------------------------------------------------------------
C     TEST FOR PROPER RANGE
C-----------------------------------------------------------------------
      AZ = ZABS(ZR,ZI)
      AA = 0.5D0/TOL
      BB=DBLE(FLOAT(I1MACH(9)))*0.5D0
      AA = DMIN1(AA,BB)
      IF (AZ.GT.AA) GO TO 260
      IF (FN.GT.AA) GO TO 260
      AA = DSQRT(AA)
      IF (AZ.GT.AA) IERR=3
      IF (FN.GT.AA) IERR=3
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
C-----------------------------------------------------------------------
      UFL = D1MACH(1)*1.0D+3
      IF (AZ.LT.UFL) GO TO 230
      IF (FNU.GT.FNUL) GO TO 90
      IF (FN.LE.1.0D0) GO TO 70
      IF (FN.GT.2.0D0) GO TO 60
      IF (AZ.GT.TOL) GO TO 70
      ARG = 0.5D0*AZ
      ALN = -FN*DLOG(ARG)
      IF (ALN.GT.ELIM) GO TO 230
      GO TO 70
   60 CONTINUE
      CALL ZUOIK(ZNR, ZNI, FNU, KODE, 2, NN, CYR, CYI, NUF, TOL, ELIM,
     * ALIM)
      IF (NUF.LT.0) GO TO 230
      NZ = NZ + NUF
      NN = NN - NUF
C-----------------------------------------------------------------------
C     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
C     IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
C-----------------------------------------------------------------------
      IF (NN.EQ.0) GO TO 140
   70 CONTINUE
      IF ((ZNR.LT.0.0D0) .OR. (ZNR.EQ.0.0D0 .AND. ZNI.LT.0.0D0 .AND.
     * M.EQ.2)) GO TO 80
C-----------------------------------------------------------------------
C     RIGHT HALF PLANE COMPUTATION, XN.GE.0. .AND. (XN.NE.0. .OR.
C     YN.GE.0. .OR. M=1)
C-----------------------------------------------------------------------
      CALL ZBKNU(ZNR, ZNI, FNU, KODE, NN, CYR, CYI, NZ, TOL, ELIM, ALIM)
      GO TO 110
C-----------------------------------------------------------------------
C     LEFT HALF PLANE COMPUTATION
C-----------------------------------------------------------------------
   80 CONTINUE
      MR = -MM
      CALL ZACON(ZNR, ZNI, FNU, KODE, MR, NN, CYR, CYI, NW, RL, FNUL,
     * TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 240
      NZ=NW
      GO TO 110
   90 CONTINUE
C-----------------------------------------------------------------------
C     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL
C-----------------------------------------------------------------------
      MR = 0
      IF ((ZNR.GE.0.0D0) .AND. (ZNR.NE.0.0D0 .OR. ZNI.GE.0.0D0 .OR.
     * M.NE.2)) GO TO 100
      MR = -MM
      IF (ZNR.NE.0.0D0 .OR. ZNI.GE.0.0D0) GO TO 100
      ZNR = -ZNR
      ZNI = -ZNI
  100 CONTINUE
      CALL ZBUNK(ZNR, ZNI, FNU, KODE, MR, NN, CYR, CYI, NW, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 240
      NZ = NZ + NW
  110 CONTINUE
C-----------------------------------------------------------------------
C     H(M,FNU,Z) = -FMM*(I/HPI)*(ZT**FNU)*K(FNU,-Z*ZT)
C
C     ZT=EXP(-FMM*HPI*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2
C-----------------------------------------------------------------------
      SGN = DSIGN(HPI,-FMM)
C-----------------------------------------------------------------------
C     CALCULATE EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      INU = INT(SNGL(FNU))
      INUH = INU/2
      IR = INU - 2*INUH
      ARG = (FNU-DBLE(FLOAT(INU-IR)))*SGN
      RHPI = 1.0D0/SGN
C     ZNI = RHPI*DCOS(ARG)
C     ZNR = -RHPI*DSIN(ARG)
      CSGNI = RHPI*DCOS(ARG)
      CSGNR = -RHPI*DSIN(ARG)
      IF (MOD(INUH,2).EQ.0) GO TO 120
C     ZNR = -ZNR
C     ZNI = -ZNI
      CSGNR = -CSGNR
      CSGNI = -CSGNI
  120 CONTINUE
      ZTI = -FMM
      RTOL = 1.0D0/TOL
      ASCLE = UFL*RTOL
      DO 130 I=1,NN
C       STR = CYR(I)*ZNR - CYI(I)*ZNI
C       CYI(I) = CYR(I)*ZNI + CYI(I)*ZNR
C       CYR(I) = STR
C       STR = -ZNI*ZTI
C       ZNI = ZNR*ZTI
C       ZNR = STR
        AA = CYR(I)
        BB = CYI(I)
        ATOL = 1.0D0
        IF (DMAX1(DABS(AA),DABS(BB)).GT.ASCLE) GO TO 135
          AA = AA*RTOL
          BB = BB*RTOL
          ATOL = TOL
  135 CONTINUE
      STR = AA*CSGNR - BB*CSGNI
      STI = AA*CSGNI + BB*CSGNR
      CYR(I) = STR*ATOL
      CYI(I) = STI*ATOL
      STR = -CSGNI*ZTI
      CSGNI = CSGNR*ZTI
      CSGNR = STR
  130 CONTINUE
      RETURN
  140 CONTINUE
      IF (ZNR.LT.0.0D0) GO TO 230
      RETURN
  230 CONTINUE
      NZ=0
      IERR=2
      RETURN
  240 CONTINUE
      IF(NW.EQ.(-1)) GO TO 230
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      NZ=0
      IERR=4
      RETURN
      END
*DECK ZBESJ
      SUBROUTINE ZBESJ (ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)
C***BEGIN PROLOGUE  ZBESJ
C***PURPOSE  Compute a sequence of the Bessel functions J(a,z) for
C            complex argument z and real nonnegative orders a=b,b+1,
C            b+2,... where b>0.  A scaling option is available to
C            help avoid overflow.
C***LIBRARY   SLATEC
C***CATEGORY  C10A4
C***TYPE      COMPLEX (CBESJ-C, ZBESJ-C)
C***KEYWORDS  BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
C             BESSEL FUNCTIONS OF THE FIRST KIND, J BESSEL FUNCTIONS
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C         On KODE=1, ZBESJ computes an N member sequence of complex
C         Bessel functions CY(L)=J(FNU+L-1,Z) for real nonnegative
C         orders FNU+L-1, L=1,...,N and complex Z in the cut plane
C         -pi<arg(Z)<=pi where Z=ZR+i*ZI.  On KODE=2, CBESJ returns
C         the scaled functions
C
C            CY(L) = exp(-abs(Y))*J(FNU+L-1,Z),  L=1,...,N and Y=Im(Z)
C
C         which remove the exponential growth in both the upper and
C         lower half planes as Z goes to infinity.  Definitions and
C         notation are found in the NBS Handbook of Mathematical
C         Functions (Ref. 1).
C
C         Input
C           ZR     - DOUBLE PRECISION real part of argument Z
C           ZI     - DOUBLE PRECISION imag part of argument Z
C           FNU    - DOUBLE PRECISION initial order, FNU>=0
C           KODE   - A parameter to indicate the scaling option
C                    KODE=1  returns
C                            CY(L)=J(FNU+L-1,Z), L=1,...,N
C                        =2  returns
C                            CY(L)=J(FNU+L-1,Z)*exp(-abs(Y)), L=1,...,N
C                            where Y=Im(Z)
C           N      - Number of terms in the sequence, N>=1
C
C         Output
C           CYR    - DOUBLE PRECISION real part of result vector
C           CYI    - DOUBLE PRECISION imag part of result vector
C           NZ     - Number of underflows set to zero
C                    NZ=0    Normal return
C                    NZ>0    CY(L)=0, L=N-NZ+1,...,N
C           IERR   - Error flag
C                    IERR=0  Normal return     - COMPUTATION COMPLETED
C                    IERR=1  Input error       - NO COMPUTATION
C                    IERR=2  Overflow          - NO COMPUTATION
C                            (Im(Z) too large on KODE=1)
C                    IERR=3  Precision warning - COMPUTATION COMPLETED
C                            (Result has half precision or less
C                            because abs(Z) or FNU+N-1 is large)
C                    IERR=4  Precision error   - NO COMPUTATION
C                            (Result has no precision because
C                            abs(Z) or FNU+N-1 is too large)
C                    IERR=5  Algorithmic error - NO COMPUTATION
C                            (Termination condition not met)
C
C *Long Description:
C
C         The computation is carried out by the formulae
C
C            J(a,z) = exp( a*pi*i/2)*I(a,-i*z),  Im(z)>=0
C
C            J(a,z) = exp(-a*pi*i/2)*I(a, i*z),  Im(z)<0
C
C         where the I Bessel function is computed as described in the
C         prologue to CBESI.
C
C         For negative orders, the formula
C
C            J(-a,z) = J(a,z)*cos(a*pi) - Y(a,z)*sin(a*pi)
C
C         can be used.  However, for large orders close to integers, the
C         the function changes radically.  When a is a large positive
C         integer, the magnitude of J(-a,z)=J(a,z)*cos(a*pi) is a
C         large negative power of ten.  But when a is not an integer,
C         Y(a,z) dominates in magnitude with a large positive power of
C         ten and the most that the second term can be reduced is by
C         unit roundoff from the coefficient.  Thus, wide changes can
C         occur within unit roundoff of a large integer for a.  Here,
C         large means a>abs(z).
C
C         In most complex variable computation, one must evaluate ele-
C         mentary functions.  When the magnitude of Z or FNU+N-1 is
C         large, losses of significance by argument reduction occur.
C         Consequently, if either one exceeds U1=SQRT(0.5/UR), then
C         losses exceeding half precision are likely and an error flag
C         IERR=3 is triggered where UR=MAX(D1MACH(4),1.0D-18) is double
C         precision unit roundoff limited to 18 digits precision.  Also,
C         if either is larger than U2=0.5/UR, then all significance is
C         lost and IERR=4.  In order to use the INT function, arguments
C         must be further restricted not to exceed the largest machine
C         integer, U3=I1MACH(9).  Thus, the magnitude of Z and FNU+N-1
C         is restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2, and
C         U3 approximate 2.0E+3, 4.2E+6, 2.1E+9 in single precision
C         and 4.7E+7, 2.3E+15 and 2.1E+9 in double precision.  This
C         makes U2 limiting in single precision and U3 limiting in
C         double precision.  This means that one can expect to retain,
C         in the worst cases on IEEE machines, no digits in single pre-
C         cision and only 6 digits in double precision.  Similar con-
C         siderations hold for other machines.
C
C         The approximate relative error in the magnitude of a complex
C         Bessel function can be expressed as P*10**S where P=MAX(UNIT
C         ROUNDOFF,1.0E-18) is the nominal precision and 10**S repre-
C         sents the increase in error due to argument reduction in the
C         elementary functions.  Here, S=MAX(1,ABS(LOG10(ABS(Z))),
C         ABS(LOG10(FNU))) approximately (i.e., S=MAX(1,ABS(EXPONENT OF
C         ABS(Z),ABS(EXPONENT OF FNU)) ).  However, the phase angle may
C         have only absolute accuracy.  This is most likely to occur
C         when one component (in magnitude) is larger than the other by
C         several orders of magnitude.  If one component is 10**K larger
C         than the other, then one can expect only MAX(ABS(LOG10(P))-K,
C         0) significant digits; or, stated another way, when K exceeds
C         the exponent of P, no significant digits remain in the smaller
C         component.  However, the phase angle retains absolute accuracy
C         because, in complex arithmetic with precision P, the smaller
C         component will not (as a rule) decrease below P times the
C         magnitude of the larger component.  In these extreme cases,
C         the principal phase angle is on the order of +P, -P, PI/2-P,
C         or -PI/2+P.
C
C***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe-
C                 matical Functions, National Bureau of Standards
C                 Applied Mathematics Series 55, U. S. Department
C                 of Commerce, Tenth Printing (1972) or later.
C               2. D. E. Amos, Computation of Bessel Functions of
C                 Complex Argument, Report SAND83-0086, Sandia National
C                 Laboratories, Albuquerque, NM, May 1983.
C               3. D. E. Amos, Computation of Bessel Functions of
C                 Complex Argument and Large Order, Report SAND83-0643,
C                 Sandia National Laboratories, Albuquerque, NM, May
C                 1983.
C               4. D. E. Amos, A Subroutine Package for Bessel Functions
C                 of a Complex Argument and Nonnegative Order, Report
C                 SAND85-1018, Sandia National Laboratory, Albuquerque,
C                 NM, May 1985.
C               5. D. E. Amos, A portable package for Bessel functions
C                 of a complex argument and nonnegative order, ACM
C                 Transactions on Mathematical Software, 12 (September
C                 1986), pp. 265-273.
C
C***ROUTINES CALLED  D1MACH, I1MACH, ZABS, ZBINU
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   890801  REVISION DATE from Version 3.2
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   920128  Category corrected.  (WRB)
C   920811  Prologue revised.  (DWL)
C***END PROLOGUE  ZBESJ
C
C     COMPLEX CI,CSGN,CY,Z,ZN
      DOUBLE PRECISION AA, ALIM, ARG, CII, CSGNI, CSGNR, CYI, CYR, DIG,
     * ELIM, FNU, FNUL, HPI, RL, R1M5, STR, TOL, ZI, ZNI, ZNR, ZR,
     * D1MACH, BB, FN, AZ, ZABS, ASCLE, RTOL, ATOL, STI
      INTEGER I, IERR, INU, INUH, IR, K, KODE, K1, K2, N, NL, NZ, I1MACH
      DIMENSION CYR(N), CYI(N)
      EXTERNAL ZABS
      DATA HPI /1.57079632679489662D0/
C
C***FIRST EXECUTABLE STATEMENT  ZBESJ
      IERR = 0
      NZ=0
      IF (FNU.LT.0.0D0) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
C-----------------------------------------------------------------------
C     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
C     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
C-----------------------------------------------------------------------
      TOL = MAX(D1MACH(4),1.0D-18)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      R1M5 = D1MACH(5)
      K = MIN(ABS(K1),ABS(K2))
      ELIM = 2.303D0*(K*R1M5-3.0D0)
      K1 = I1MACH(14) - 1
      AA = R1M5*K1
      DIG = MIN(AA,18.0D0)
      AA = AA*2.303D0
      ALIM = ELIM + MAX(-AA,-41.45D0)
      RL = 1.2D0*DIG + 3.0D0
      FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
C-----------------------------------------------------------------------
C     TEST FOR PROPER RANGE
C-----------------------------------------------------------------------
      AZ = ZABS(ZR,ZI)
      FN = FNU+(N-1)
      AA = 0.5D0/TOL
      BB = I1MACH(9)*0.5D0
      AA = MIN(AA,BB)
      IF (AZ.GT.AA) GO TO 260
      IF (FN.GT.AA) GO TO 260
      AA = SQRT(AA)
      IF (AZ.GT.AA) IERR=3
      IF (FN.GT.AA) IERR=3
C-----------------------------------------------------------------------
C     CALCULATE CSGN=EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C     WHEN FNU IS LARGE
C-----------------------------------------------------------------------
      CII = 1.0D0
      INU = FNU
      INUH = INU/2
      IR = INU - 2*INUH
      ARG = (FNU-(INU-IR))*HPI
      CSGNR = COS(ARG)
      CSGNI = SIN(ARG)
      IF (MOD(INUH,2).EQ.0) GO TO 40
      CSGNR = -CSGNR
      CSGNI = -CSGNI
   40 CONTINUE
C-----------------------------------------------------------------------
C     ZN IS IN THE RIGHT HALF PLANE
C-----------------------------------------------------------------------
      ZNR = ZI
      ZNI = -ZR
      IF (ZI.GE.0.0D0) GO TO 50
      ZNR = -ZNR
      ZNI = -ZNI
      CSGNI = -CSGNI
      CII = -CII
   50 CONTINUE
      CALL ZBINU(ZNR, ZNI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL, TOL,
     * ELIM, ALIM)
      IF (NZ.LT.0) GO TO 130
      NL = N - NZ
      IF (NL.EQ.0) RETURN
      RTOL = 1.0D0/TOL
      ASCLE = D1MACH(1)*RTOL*1.0D+3
      DO 60 I=1,NL
C       STR = CYR(I)*CSGNR - CYI(I)*CSGNI
C       CYI(I) = CYR(I)*CSGNI + CYI(I)*CSGNR
C       CYR(I) = STR
        AA = CYR(I)
        BB = CYI(I)
        ATOL = 1.0D0
        IF (MAX(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 55
          AA = AA*RTOL
          BB = BB*RTOL
          ATOL = TOL
   55   CONTINUE
        STR = AA*CSGNR - BB*CSGNI
        STI = AA*CSGNI + BB*CSGNR
        CYR(I) = STR*ATOL
        CYI(I) = STI*ATOL
        STR = -CSGNI*CII
        CSGNI = CSGNR*CII
        CSGNR = STR
   60 CONTINUE
      RETURN
  130 CONTINUE
      IF(NZ.EQ.(-2)) GO TO 140
      NZ = 0
      IERR = 2
      RETURN
  140 CONTINUE
      NZ=0
      IERR=5
      RETURN
  260 CONTINUE
      NZ=0
      IERR=4
      RETURN
      END
*DECK ZBESY
      SUBROUTINE ZBESY (ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, CWRKR,
     +   CWRKI, IERR)
C***BEGIN PROLOGUE  ZBESY
C***PURPOSE  Compute a sequence of the Bessel functions Y(a,z) for
C            complex argument z and real nonnegative orders a=b,b+1,
C            b+2,... where b>0.  A scaling option is available to
C            help avoid overflow.
C***LIBRARY   SLATEC
C***CATEGORY  C10A4
C***TYPE      COMPLEX (CBESY-C, ZBESY-C)
C***KEYWORDS  BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
C             BESSEL FUNCTIONS OF SECOND KIND, WEBER'S FUNCTION,
C             Y BESSEL FUNCTIONS
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C                      ***A DOUBLE PRECISION ROUTINE***
C         On KODE=1, ZBESY computes an N member sequence of complex
C         Bessel functions CY(L)=Y(FNU+L-1,Z) for real nonnegative
C         orders FNU+L-1, L=1,...,N and complex Z in the cut plane
C         -pi<arg(Z)<=pi where Z=ZR+i*ZI.  On KODE=2, CBESY returns
C         the scaled functions
C
C            CY(L) = exp(-abs(Y))*Y(FNU+L-1,Z),  L=1,...,N, Y=Im(Z)
C
C         which remove the exponential growth in both the upper and
C         lower half planes as Z goes to infinity.  Definitions and
C         notation are found in the NBS Handbook of Mathematical
C         Functions (Ref. 1).
C
C         Input
C           ZR     - DOUBLE PRECISION real part of nonzero argument Z
C           ZI     - DOUBLE PRECISION imag part of nonzero argument Z
C           FNU    - DOUBLE PRECISION initial order, FNU>=0
C           KODE   - A parameter to indicate the scaling option
C                    KODE=1  returns
C                            CY(L)=Y(FNU+L-1,Z), L=1,...,N
C                        =2  returns
C                            CY(L)=Y(FNU+L-1,Z)*exp(-abs(Y)), L=1,...,N
C                            where Y=Im(Z)
C           N      - Number of terms in the sequence, N>=1
C           CWRKR  - DOUBLE PRECISION work vector of dimension N
C           CWRKI  - DOUBLE PRECISION work vector of dimension N
C
C         Output
C           CYR    - DOUBLE PRECISION real part of result vector
C           CYI    - DOUBLE PRECISION imag part of result vector
C           NZ     - Number of underflows set to zero
C                    NZ=0    Normal return
C                    NZ>0    CY(L)=0 for NZ values of L, usually on
C                            KODE=2 (the underflows may not be in an
C                            uninterrupted sequence)
C           IERR   - Error flag
C                    IERR=0  Normal return     - COMPUTATION COMPLETED
C                    IERR=1  Input error       - NO COMPUTATION
C                    IERR=2  Overflow          - NO COMPUTATION
C                            (abs(Z) too small and/or FNU+N-1
C                            too large)
C                    IERR=3  Precision warning - COMPUTATION COMPLETED
C                            (Result has half precision or less
C                            because abs(Z) or FNU+N-1 is large)
C                    IERR=4  Precision error   - NO COMPUTATION
C                            (Result has no precision because
C                            abs(Z) or FNU+N-1 is too large)
C                    IERR=5  Algorithmic error - NO COMPUTATION
C                            (Termination condition not met)
C
C *Long Description:
C
C         The computation is carried out by the formula
C
C            Y(a,z) = (H(1,a,z) - H(2,a,z))/(2*i)
C
C         where the Hankel functions are computed as described in CBESH.
C
C         For negative orders, the formula
C
C            Y(-a,z) = Y(a,z)*cos(a*pi) + J(a,z)*sin(a*pi)
C
C         can be used.  However, for large orders close to half odd
C         integers the function changes radically.  When a is a large
C         positive half odd integer, the magnitude of Y(-a,z)=J(a,z)*
C         sin(a*pi) is a large negative power of ten.  But when a is
C         not a half odd integer, Y(a,z) dominates in magnitude with a
C         large positive power of ten and the most that the second term
C         can be reduced is by unit roundoff from the coefficient.
C         Thus,  wide changes can occur within unit roundoff of a large
C         half odd integer.  Here, large means a>abs(z).
C
C         In most complex variable computation, one must evaluate ele-
C         mentary functions.  When the magnitude of Z or FNU+N-1 is
C         large, losses of significance by argument reduction occur.
C         Consequently, if either one exceeds U1=SQRT(0.5/UR), then
C         losses exceeding half precision are likely and an error flag
C         IERR=3 is triggered where UR=MAX(D1MACH(4),1.0D-18) is double
C         precision unit roundoff limited to 18 digits precision.  Also,
C         if either is larger than U2=0.5/UR, then all significance is
C         lost and IERR=4.  In order to use the INT function, arguments
C         must be further restricted not to exceed the largest machine
C         integer, U3=I1MACH(9).  Thus, the magnitude of Z and FNU+N-1
C         is restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2, and
C         U3 approximate 2.0E+3, 4.2E+6, 2.1E+9 in single precision
C         and 4.7E+7, 2.3E+15 and 2.1E+9 in double precision.  This
C         makes U2 limiting in single precision and U3 limiting in
C         double precision.  This means that one can expect to retain,
C         in the worst cases on IEEE machines, no digits in single pre-
C         cision and only 6 digits in double precision.  Similar con-
C         siderations hold for other machines.
C
C         The approximate relative error in the magnitude of a complex
C         Bessel function can be expressed as P*10**S where P=MAX(UNIT
C         ROUNDOFF,1.0E-18) is the nominal precision and 10**S repre-
C         sents the increase in error due to argument reduction in the
C         elementary functions.  Here, S=MAX(1,ABS(LOG10(ABS(Z))),
C         ABS(LOG10(FNU))) approximately (i.e., S=MAX(1,ABS(EXPONENT OF
C         ABS(Z),ABS(EXPONENT OF FNU)) ).  However, the phase angle may
C         have only absolute accuracy.  This is most likely to occur
C         when one component (in magnitude) is larger than the other by
C         several orders of magnitude.  If one component is 10**K larger
C         than the other, then one can expect only MAX(ABS(LOG10(P))-K,
C         0) significant digits; or, stated another way, when K exceeds
C         the exponent of P, no significant digits remain in the smaller
C         component.  However, the phase angle retains absolute accuracy
C         because, in complex arithmetic with precision P, the smaller
C         component will not (as a rule) decrease below P times the
C         magnitude of the larger component.  In these extreme cases,
C         the principal phase angle is on the order of +P, -P, PI/2-P,
C         or -PI/2+P.
C
C***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe-
C                 matical Functions, National Bureau of Standards
C                 Applied Mathematics Series 55, U. S. Department
C                 of Commerce, Tenth Printing (1972) or later.
C               2. D. E. Amos, Computation of Bessel Functions of
C                 Complex Argument, Report SAND83-0086, Sandia National
C                 Laboratories, Albuquerque, NM, May 1983.
C               3. D. E. Amos, Computation of Bessel Functions of
C                 Complex Argument and Large Order, Report SAND83-0643,
C                 Sandia National Laboratories, Albuquerque, NM, May
C                 1983.
C               4. D. E. Amos, A Subroutine Package for Bessel Functions
C                 of a Complex Argument and Nonnegative Order, Report
C                 SAND85-1018, Sandia National Laboratory, Albuquerque,
C                 NM, May 1985.
C               5. D. E. Amos, A portable package for Bessel functions
C                 of a complex argument and nonnegative order, ACM
C                 Transactions on Mathematical Software, 12 (September
C                 1986), pp. 265-273.
C
C***ROUTINES CALLED  D1MACH, I1MACH, ZBESH
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   890801  REVISION DATE from Version 3.2
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   920128  Category corrected.  (WRB)
C   920811  Prologue revised.  (DWL)
C***END PROLOGUE  ZBESY
C
C     COMPLEX CWRK,CY,C1,C2,EX,HCI,Z,ZU,ZV
      DOUBLE PRECISION CWRKI, CWRKR, CYI, CYR, C1I, C1R, C2I, C2R,
     * ELIM, EXI, EXR, EY, FNU, HCII, STI, STR, TAY, ZI, ZR,
     * D1MACH, ASCLE, RTOL, ATOL, AA, BB, TOL, R1M5
      INTEGER I, IERR, K, KODE, K1, K2, N, NZ, NZ1, NZ2, I1MACH
      DIMENSION CYR(N), CYI(N), CWRKR(N), CWRKI(N)
C***FIRST EXECUTABLE STATEMENT  ZBESY
      IERR = 0
      NZ=0
      IF (ZR.EQ.0.0D0 .AND. ZI.EQ.0.0D0) IERR=1
      IF (FNU.LT.0.0D0) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (N.LT.1) IERR=1
      IF (IERR.NE.0) RETURN
      HCII = 0.5D0
      CALL ZBESH(ZR, ZI, FNU, KODE, 1, N, CYR, CYI, NZ1, IERR)
      IF (IERR.NE.0.AND.IERR.NE.3) GO TO 170
      CALL ZBESH(ZR, ZI, FNU, KODE, 2, N, CWRKR, CWRKI, NZ2, IERR)
      IF (IERR.NE.0.AND.IERR.NE.3) GO TO 170
      NZ = MIN(NZ1,NZ2)
      IF (KODE.EQ.2) GO TO 60
      DO 50 I=1,N
        STR = CWRKR(I) - CYR(I)
        STI = CWRKI(I) - CYI(I)
        CYR(I) = -STI*HCII
        CYI(I) = STR*HCII
   50 CONTINUE
      RETURN
   60 CONTINUE
      TOL = MAX(D1MACH(4),1.0D-18)
      K1 = I1MACH(15)
      K2 = I1MACH(16)
      K = MIN(ABS(K1),ABS(K2))
      R1M5 = D1MACH(5)
C-----------------------------------------------------------------------
C     ELIM IS THE APPROXIMATE EXPONENTIAL UNDER- AND OVERFLOW LIMIT
C-----------------------------------------------------------------------
      ELIM = 2.303D0*(K*R1M5-3.0D0)
      EXR = COS(ZR)
      EXI = SIN(ZR)
      EY = 0.0D0
      TAY = ABS(ZI+ZI)
      IF (TAY.LT.ELIM) EY = EXP(-TAY)
      IF (ZI.LT.0.0D0) GO TO 90
      C1R = EXR*EY
      C1I = EXI*EY
      C2R = EXR
      C2I = -EXI
   70 CONTINUE
      NZ = 0
      RTOL = 1.0D0/TOL
      ASCLE = D1MACH(1)*RTOL*1.0D+3
      DO 80 I=1,N
C       STR = C1R*CYR(I) - C1I*CYI(I)
C       STI = C1R*CYI(I) + C1I*CYR(I)
C       STR = -STR + C2R*CWRKR(I) - C2I*CWRKI(I)
C       STI = -STI + C2R*CWRKI(I) + C2I*CWRKR(I)
C       CYR(I) = -STI*HCII
C       CYI(I) = STR*HCII
        AA = CWRKR(I)
        BB = CWRKI(I)
        ATOL = 1.0D0
        IF (MAX(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 75
          AA = AA*RTOL
          BB = BB*RTOL
          ATOL = TOL
   75   CONTINUE
        STR = (AA*C2R - BB*C2I)*ATOL
        STI = (AA*C2I + BB*C2R)*ATOL
        AA = CYR(I)
        BB = CYI(I)
        ATOL = 1.0D0
        IF (MAX(ABS(AA),ABS(BB)).GT.ASCLE) GO TO 85
          AA = AA*RTOL
          BB = BB*RTOL
          ATOL = TOL
   85   CONTINUE
        STR = STR - (AA*C1R - BB*C1I)*ATOL
        STI = STI - (AA*C1I + BB*C1R)*ATOL
        CYR(I) = -STI*HCII
        CYI(I) =  STR*HCII
        IF (STR.EQ.0.0D0 .AND. STI.EQ.0.0D0 .AND. EY.EQ.0.0D0) NZ = NZ
     *   + 1
   80 CONTINUE
      RETURN
   90 CONTINUE
      C1R = EXR
      C1I = EXI
      C2R = EXR*EY
      C2I = -EXI*EY
      GO TO 70
  170 CONTINUE
      NZ = 0
      RETURN
      END
*DECK ZBINU
      SUBROUTINE ZBINU (ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, RL, FNUL,
     +   TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  ZBINU
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK and ZBIRY
C***LIBRARY   SLATEC
C***TYPE      ALL (CBINU-A, ZBINU-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE
C
C***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBIRY
C***ROUTINES CALLED  ZABS, ZASYI, ZBUNI, ZMLRI, ZSERI, ZUOIK, ZWRSK
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZBINU
      DOUBLE PRECISION ALIM, AZ, CWI, CWR, CYI, CYR, DFNU, ELIM, FNU,
     * FNUL, RL, TOL, ZEROI, ZEROR, ZI, ZR, ZABS
      INTEGER I, INW, KODE, N, NLAST, NN, NUI, NW, NZ
      DIMENSION CYR(N), CYI(N), CWR(2), CWI(2)
      EXTERNAL ZABS
      DATA ZEROR,ZEROI / 0.0D0, 0.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZBINU
      NZ = 0
      AZ = ZABS(ZR,ZI)
      NN = N
      DFNU = FNU + (N-1)
      IF (AZ.LE.2.0D0) GO TO 10
      IF (AZ*AZ*0.25D0.GT.DFNU+1.0D0) GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     POWER SERIES
C-----------------------------------------------------------------------
      CALL ZSERI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL, ELIM, ALIM)
      INW = ABS(NW)
      NZ = NZ + INW
      NN = NN - INW
      IF (NN.EQ.0) RETURN
      IF (NW.GE.0) GO TO 120
      DFNU = FNU + (NN-1)
   20 CONTINUE
      IF (AZ.LT.RL) GO TO 40
      IF (DFNU.LE.1.0D0) GO TO 30
      IF (AZ+AZ.LT.DFNU*DFNU) GO TO 50
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR LARGE Z
C-----------------------------------------------------------------------
   30 CONTINUE
      CALL ZASYI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, RL, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 130
      GO TO 120
   40 CONTINUE
      IF (DFNU.LE.1.0D0) GO TO 70
   50 CONTINUE
C-----------------------------------------------------------------------
C     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
C-----------------------------------------------------------------------
      CALL ZUOIK(ZR, ZI, FNU, KODE, 1, NN, CYR, CYI, NW, TOL, ELIM,
     * ALIM)
      IF (NW.LT.0) GO TO 130
      NZ = NZ + NW
      NN = NN - NW
      IF (NN.EQ.0) RETURN
      DFNU = FNU+(NN-1)
      IF (DFNU.GT.FNUL) GO TO 110
      IF (AZ.GT.FNUL) GO TO 110
   60 CONTINUE
      IF (AZ.GT.RL) GO TO 80
   70 CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE SERIES
C-----------------------------------------------------------------------
      CALL ZMLRI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, TOL)
      IF(NW.LT.0) GO TO 130
      GO TO 120
   80 CONTINUE
C-----------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
C-----------------------------------------------------------------------
      CALL ZUOIK(ZR, ZI, FNU, KODE, 2, 2, CWR, CWI, NW, TOL, ELIM,
     * ALIM)
      IF (NW.GE.0) GO TO 100
      NZ = NN
      DO 90 I=1,NN
        CYR(I) = ZEROR
        CYI(I) = ZEROI
   90 CONTINUE
      RETURN
  100 CONTINUE
      IF (NW.GT.0) GO TO 130
      CALL ZWRSK(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, CWR, CWI, TOL,
     * ELIM, ALIM)
      IF (NW.LT.0) GO TO 130
      GO TO 120
  110 CONTINUE
C-----------------------------------------------------------------------
C     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
C-----------------------------------------------------------------------
      NUI = FNUL-DFNU + 1
      NUI = MAX(NUI,0)
      CALL ZBUNI(ZR, ZI, FNU, KODE, NN, CYR, CYI, NW, NUI, NLAST, FNUL,
     * TOL, ELIM, ALIM)
      IF (NW.LT.0) GO TO 130
      NZ = NZ + NW
      IF (NLAST.EQ.0) GO TO 120
      NN = NLAST
      GO TO 60
  120 CONTINUE
      RETURN
  130 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END
      subroutine zbjy(x,m,nstop,fnu,
     & 							bjr,byr,cjr,cji)
c
c *********************************************************************
c
c  Subroutine zbjy gets J & Y Bessel functions for use in
c  sphere (fnu=0.50) or cylinder (fnu=0.0d0) scattering calculations.
c  Scaled or nonscaled functions for argument z = m*x are returned
c  depending
c  on the magnitude of the product of the size parameter and the complex
c  refractive index (zabs(z)). Nonscaled functions are returned if the
c  imaginary part of the refractive index is zero. Nonscaled functions are
c  returned for argument x.
c
c       Merrill Milham      >>> version: 2.0 <<<        JANUARY 1994
C
c  inputs:
c       x = the size parameter of the cylinder.(real*8)
c       m = the complex refractive index, n - ik.(complex*16)
c   nstop = the highest order of the Bessel functions.(integer)
c     fnu = 0.5d0 for sphere calculations or
c           0.0d0 for cylinder calculations. (real*8)
c
c  outputs:
c
c   bjr = real part of J(x) Bessel functions.(array: real*8)
c   byr = real part of Y(x) Bessel functions.(array: real*8)
c   cjr = real part of J(m*x) Bessel functions.(array: real*8)
c   cji = imag. part of J(m*x) Bessel functions.(array: real*8)
c
c subroutines used:
c
c           zbesj - returns J Bessel functions
c           zbesy - returns Y Bessel functions
c
c Reference: D. E. Amos, "Algorthim 644: A Portable Package for Bessel
c            Functions of a Complex Argument and Nonnegative Order,"
c            ACM Transcations on Mathematical Software, 12,265-273(1986).
c
c *************************************************************************
c
      implicit none
cthw      automatic cwrkr,cwrki,bji,byi
*
      integer al
      real*8 zero,xll,two
      parameter (al=5100,zero=0.0d0,xll=1.d-1,two=2.d0) 
*
      real*8 bjr(1) ,byr(1)
c
c page 31 
c
      real*8 x,fnu,cjr(1),cji(1)
      real*8 zr,zi,d1mach
      real*8 cwrkr(al),cwrki(al)
      real*8 r1m5,elim,aa,alim
      complex*16 m,z
      integer kode,ierr,nz,nstop,n
      integer k1,i1mach,k2,k
      logical zflag,eflg1,eflg2
c
      kode=1
c
      call zbesj(x,zero,fnu,kode,nstop,bjr,cji,nz,ierr) 
*
      eflg1=ierr.eq.0
      eflg2=nz.eq.0
      if (eflg1.and.eflg2) then
            continue
        else if (.not.eflg2.and.eflg1) then
                  nstop=nstop-nz
                  print*,'zbesj error: ierr =',ierr,'nz =',nz
                    else
              print*,'inputs =',x,zero,fnu,kode,nstop
              print*,'zbesj called from subroutine zbjy'
      end if
c
      call zbesy(x,zero,fnu,kode,nstop,byr,cji,nz,cwrkr,cwrki,ierr) 
*
      eflg1=ierr.eq.0
      eflg2=nz.eq.0
      if (eflg1.and.eflg2) then
            continue
            else if (.not.eflg2.and.eflg1) then
                        nstop=nstop-nz
                        print*,'zbesy error: ierr =',ierr,'nz =',nz
            continue
                    else
                print*,'zbesy error: ierr =',ierr,'nz=',nz
                print*,'inputs =',x,zero,fnu,kode,nstop
                print*,'zbesy called from subroutine zbjy'
      end if
C
      z=m*x
      zr=dble(z)
      zi=dimag(z)
      if(zi.ne.zero) then 
*
      k1 = i1mach(15)
      k2 = i1mach(16)
      r1m5 = d1mach(5)
      k = min0(iabs(k1),iabs(k2))
      elim = 2.303d0*(dble(float(k))*r1m5-3.0d0)
      k1 = i1mach(14) - 1
c
c page 32 
c
      aa = r1m5*dble(float(k1))
      aa = aa*2.303d0
      alim = elim + dmax1(-aa,-41.45d0) 
*
            if(dabs(zi).gt.alim) kode=2
                                    else
           continue
      end if
      call zbesj(zr,zi,fnu,kode,nstop,cjr,cji,nz,ierr)
*
      eflg1=ierr.eq.0
      eflg2=nz.eq.0
      if (eflg1.and.eflg2) then
            continue
            else if (.not.eflg2.and.eflg1) then
                        nstop=nstop-nz
                        print*,'zbesj error: ierr =',ierr,'nz =',nz
            continue
                        else
                print*,'zbesj error: ierr =',ierr,'nz =',nz
                print*,'inputs =',zr,zi,fnu,kode,nstop
                print*,'zbesj called from subroutine zbjy'
      end if
c
      zflag=dabs(zi).lt.two*d1mach(1).and.x.lt.xll
      if (.not.zflag) then
            continue
                    else
            do 100 n=1,nstop
                cji (n)=zero
100         continue
      end if 
*
      zflag=dabs(zr).lt.two*d1mach(1).and.x.lt.xll
      if (.not.zflag) then
            continue
                      else
            cji(1)=zero
            do 200 n=2,nstop
                    if(mod(n,2).eq.0) then
                        cjr(n)=zero
                                       else
                        cji(n)=zero
                    end if
200     continue
      end if
c
      return
c
      end
c*DECK ZBKNU
      SUBROUTINE ZBKNU (ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM,
     +   ALIM)
C***BEGIN PROLOGUE  ZBKNU
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZAIRY, ZBESH, ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CBKNU-A, ZBKNU-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE.
C
C***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, DGAMLN, I1MACH, ZABS, ZDIV, ZEXP, ZKSCL,
C                    ZLOG, ZMLT, ZSHCH, ZSQRT, ZUCHK
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   930122  Added ZEXP, ZLOG and ZSQRT to EXTERNAL statement.  (RWC)
C***END PROLOGUE  ZBKNU
C
      DOUBLE PRECISION AA, AK, ALIM, ASCLE, A1, A2, BB, BK, BRY, CAZ,
     * CBI, CBR, CC, CCHI, CCHR, CKI, CKR, COEFI, COEFR, CONEI, CONER,
     * CRSCR, CSCLR, CSHI, CSHR, CSI, CSR, CSRR, CSSR, CTWOR,
     * CZEROI, CZEROR, CZI, CZR, DNU, DNU2, DPI, ELIM, ETEST, FC, FHS,
     * FI, FK, FKS, FMUI, FMUR, FNU, FPI, FR, G1, G2, HPI, PI, PR, PTI,
     * PTR, P1I, P1R, P2I, P2M, P2R, QI, QR, RAK, RCAZ, RTHPI, RZI,
     * RZR, R1, S, SMUI, SMUR, SPI, STI, STR, S1I, S1R, S2I, S2R, TM,
     * TOL, TTH, T1, T2, YI, YR, ZI, ZR, DGAMLN, D1MACH, ZABS, ELM,
     * CELMR, ZDR, ZDI, AS, ALAS, HELIM, CYR, CYI
      INTEGER I, IFLAG, INU, K, KFLAG, KK, KMAX, KODE, KODED, N, NZ,
     * IDUM, I1MACH, J, IC, INUB, NW
      DIMENSION YR(N), YI(N), CC(8), CSSR(3), CSRR(3), BRY(3), CYR(2),
     * CYI(2)
      EXTERNAL ZABS, ZEXP, ZLOG, ZSQRT
C     COMPLEX Z,Y,A,B,RZ,SMU,FU,FMU,F,FLRZ,CZ,S1,S2,CSH,CCH
C     COMPLEX CK,P,Q,COEF,P1,P2,CBK,PT,CZERO,CONE,CTWO,ST,EZ,CS,DK
C
      DATA KMAX / 30 /
      DATA CZEROR,CZEROI,CONER,CONEI,CTWOR,R1/
     1  0.0D0 , 0.0D0 , 1.0D0 , 0.0D0 , 2.0D0 , 2.0D0 /
      DATA DPI, RTHPI, SPI ,HPI, FPI, TTH /
     1     3.14159265358979324D0,       1.25331413731550025D0,
     2     1.90985931710274403D0,       1.57079632679489662D0,
     3     1.89769999331517738D0,       6.66666666666666666D-01/
      DATA CC(1), CC(2), CC(3), CC(4), CC(5), CC(6), CC(7), CC(8)/
     1     5.77215664901532861D-01,    -4.20026350340952355D-02,
     2    -4.21977345555443367D-02,     7.21894324666309954D-03,
     3    -2.15241674114950973D-04,    -2.01348547807882387D-05,
     4     1.13302723198169588D-06,     6.11609510448141582D-09/
C***FIRST EXECUTABLE STATEMENT  ZBKNU
      CAZ = ZABS(ZR,ZI)
      CSCLR = 1.0D0/TOL
      CRSCR = TOL
      CSSR(1) = CSCLR
      CSSR(2) = 1.0D0
      CSSR(3) = CRSCR
      CSRR(1) = CRSCR
      CSRR(2) = 1.0D0
      CSRR(3) = CSCLR
      BRY(1) = 1.0D+3*D1MACH(1)/TOL
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = D1MACH(2)
      NZ = 0
      IFLAG = 0
      KODED = KODE
      RCAZ = 1.0D0/CAZ
      STR = ZR*RCAZ
      STI = -ZI*RCAZ
      RZR = (STR+STR)*RCAZ
      RZI = (STI+STI)*RCAZ
      INU = FNU+0.5D0
      DNU = FNU - INU
      IF (ABS(DNU).EQ.0.5D0) GO TO 110
      DNU2 = 0.0D0
      IF (ABS(DNU).GT.TOL) DNU2 = DNU*DNU
      IF (CAZ.GT.R1) GO TO 110
C-----------------------------------------------------------------------
C     SERIES FOR ABS(Z).LE.R1
C-----------------------------------------------------------------------
      FC = 1.0D0
      CALL ZLOG(RZR, RZI, SMUR, SMUI, IDUM)
      FMUR = SMUR*DNU
      FMUI = SMUI*DNU
      CALL ZSHCH(FMUR, FMUI, CSHR, CSHI, CCHR, CCHI)
      IF (DNU.EQ.0.0D0) GO TO 10
      FC = DNU*DPI
      FC = FC/SIN(FC)
      SMUR = CSHR/DNU
      SMUI = CSHI/DNU
   10 CONTINUE
      A2 = 1.0D0 + DNU
C-----------------------------------------------------------------------
C     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU)
C-----------------------------------------------------------------------
      T2 = EXP(-DGAMLN(A2,IDUM))
      T1 = 1.0D0/(T2*FC)
      IF (ABS(DNU).GT.0.1D0) GO TO 40
C-----------------------------------------------------------------------
C     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
C-----------------------------------------------------------------------
      AK = 1.0D0
      S = CC(1)
      DO 20 K=2,8
        AK = AK*DNU2
        TM = CC(K)*AK
        S = S + TM
        IF (ABS(TM).LT.TOL) GO TO 30
   20 CONTINUE
   30 G1 = -S
      GO TO 50
   40 CONTINUE
      G1 = (T1-T2)/(DNU+DNU)
   50 CONTINUE
      G2 = (T1+T2)*0.5D0
      FR = FC*(CCHR*G1+SMUR*G2)
      FI = FC*(CCHI*G1+SMUI*G2)
      CALL ZEXP(FMUR, FMUI, STR, STI)
      PR = 0.5D0*STR/T2
      PI = 0.5D0*STI/T2
      CALL ZDIV(0.5D0, 0.0D0, STR, STI, PTR, PTI)
      QR = PTR/T1
      QI = PTI/T1
      S1R = FR
      S1I = FI
      S2R = PR
      S2I = PI
      AK = 1.0D0
      A1 = 1.0D0
      CKR = CONER
      CKI = CONEI
      BK = 1.0D0 - DNU2
      IF (INU.GT.0 .OR. N.GT.1) GO TO 80
C-----------------------------------------------------------------------
C     GENERATE K(FNU,Z), 0.0D0 .LE. FNU .LT. 0.5D0 AND N=1
C-----------------------------------------------------------------------
      IF (CAZ.LT.TOL) GO TO 70
      CALL ZMLT(ZR, ZI, ZR, ZI, CZR, CZI)
      CZR = 0.25D0*CZR
      CZI = 0.25D0*CZI
      T1 = 0.25D0*CAZ*CAZ
   60 CONTINUE
      FR = (FR*AK+PR+QR)/BK
      FI = (FI*AK+PI+QI)/BK
      STR = 1.0D0/(AK-DNU)
      PR = PR*STR
      PI = PI*STR
      STR = 1.0D0/(AK+DNU)
      QR = QR*STR
      QI = QI*STR
      STR = CKR*CZR - CKI*CZI
      RAK = 1.0D0/AK
      CKI = (CKR*CZI+CKI*CZR)*RAK
      CKR = STR*RAK
      S1R = CKR*FR - CKI*FI + S1R
      S1I = CKR*FI + CKI*FR + S1I
      A1 = A1*T1*RAK
      BK = BK + AK + AK + 1.0D0
      AK = AK + 1.0D0
      IF (A1.GT.TOL) GO TO 60
   70 CONTINUE
      YR(1) = S1R
      YI(1) = S1I
      IF (KODED.EQ.1) RETURN
      CALL ZEXP(ZR, ZI, STR, STI)
      CALL ZMLT(S1R, S1I, STR, STI, YR(1), YI(1))
      RETURN
C-----------------------------------------------------------------------
C     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
C-----------------------------------------------------------------------
   80 CONTINUE
      IF (CAZ.LT.TOL) GO TO 100
      CALL ZMLT(ZR, ZI, ZR, ZI, CZR, CZI)
      CZR = 0.25D0*CZR
      CZI = 0.25D0*CZI
      T1 = 0.25D0*CAZ*CAZ
   90 CONTINUE
      FR = (FR*AK+PR+QR)/BK
      FI = (FI*AK+PI+QI)/BK
      STR = 1.0D0/(AK-DNU)
      PR = PR*STR
      PI = PI*STR
      STR = 1.0D0/(AK+DNU)
      QR = QR*STR
      QI = QI*STR
      STR = CKR*CZR - CKI*CZI
      RAK = 1.0D0/AK
      CKI = (CKR*CZI+CKI*CZR)*RAK
      CKR = STR*RAK
      S1R = CKR*FR - CKI*FI + S1R
      S1I = CKR*FI + CKI*FR + S1I
      STR = PR - FR*AK
      STI = PI - FI*AK
      S2R = CKR*STR - CKI*STI + S2R
      S2I = CKR*STI + CKI*STR + S2I
      A1 = A1*T1*RAK
      BK = BK + AK + AK + 1.0D0
      AK = AK + 1.0D0
      IF (A1.GT.TOL) GO TO 90
  100 CONTINUE
      KFLAG = 2
      A1 = FNU + 1.0D0
      AK = A1*ABS(SMUR)
      IF (AK.GT.ALIM) KFLAG = 3
      STR = CSSR(KFLAG)
      P2R = S2R*STR
      P2I = S2I*STR
      CALL ZMLT(P2R, P2I, RZR, RZI, S2R, S2I)
      S1R = S1R*STR
      S1I = S1I*STR
      IF (KODED.EQ.1) GO TO 210
      CALL ZEXP(ZR, ZI, FR, FI)
      CALL ZMLT(S1R, S1I, FR, FI, S1R, S1I)
      CALL ZMLT(S2R, S2I, FR, FI, S2R, S2I)
      GO TO 210
C-----------------------------------------------------------------------
C     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
C     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
C     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
C     RECURSION
C-----------------------------------------------------------------------
  110 CONTINUE
      CALL ZSQRT(ZR, ZI, STR, STI)
      CALL ZDIV(RTHPI, CZEROI, STR, STI, COEFR, COEFI)
      KFLAG = 2
      IF (KODED.EQ.2) GO TO 120
      IF (ZR.GT.ALIM) GO TO 290
C     BLANK LINE
      STR = EXP(-ZR)*CSSR(KFLAG)
      STI = -STR*SIN(ZI)
      STR = STR*COS(ZI)
      CALL ZMLT(COEFR, COEFI, STR, STI, COEFR, COEFI)
  120 CONTINUE
      IF (ABS(DNU).EQ.0.5D0) GO TO 300
C-----------------------------------------------------------------------
C     MILLER ALGORITHM FOR ABS(Z).GT.R1
C-----------------------------------------------------------------------
      AK = COS(DPI*DNU)
      AK = ABS(AK)
      IF (AK.EQ.CZEROR) GO TO 300
      FHS = ABS(0.25D0-DNU2)
      IF (FHS.EQ.CZEROR) GO TO 300
C-----------------------------------------------------------------------
C     COMPUTE R2=F(E). IF ABS(Z).GE.R2, USE FORWARD RECURRENCE TO
C     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
C     12.LE.E.LE.60. E IS COMPUTED FROM 2**(-E)=B**(1-I1MACH(14))=
C     TOL WHERE B IS THE BASE OF THE ARITHMETIC.
C-----------------------------------------------------------------------
      T1 = I1MACH(14)-1
      T1 = T1*D1MACH(5)*3.321928094D0
      T1 = MAX(T1,12.0D0)
      T1 = MIN(T1,60.0D0)
      T2 = TTH*T1 - 6.0D0
      IF (ZR.NE.0.0D0) GO TO 130
      T1 = HPI
      GO TO 140
  130 CONTINUE
      T1 = DATAN(ZI/ZR)
      T1 = ABS(T1)
  140 CONTINUE
      IF (T2.GT.CAZ) GO TO 170
C-----------------------------------------------------------------------
C     FORWARD RECURRENCE LOOP WHEN ABS(Z).GE.R2
C-----------------------------------------------------------------------
      ETEST = AK/(DPI*CAZ*TOL)
      FK = CONER
      IF (ETEST.LT.CONER) GO TO 180
      FKS = CTWOR
      CKR = CAZ + CAZ + CTWOR
      P1R = CZEROR
      P2R = CONER
      DO 150 I=1,KMAX
        AK = FHS/FKS
        CBR = CKR/(FK+CONER)
        PTR = P2R
        P2R = CBR*P2R - P1R*AK
        P1R = PTR
        CKR = CKR + CTWOR
        FKS = FKS + FK + FK + CTWOR
        FHS = FHS + FK + FK
        FK = FK + CONER
        STR = ABS(P2R)*FK
        IF (ETEST.LT.STR) GO TO 160
  150 CONTINUE
      GO TO 310
  160 CONTINUE
      FK = FK + SPI*T1*SQRT(T2/CAZ)
      FHS = ABS(0.25D0-DNU2)
      GO TO 180
  170 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE BACKWARD INDEX K FOR ABS(Z).LT.R2
C-----------------------------------------------------------------------
      A2 = SQRT(CAZ)
      AK = FPI*AK/(TOL*SQRT(A2))
      AA = 3.0D0*T1/(1.0D0+CAZ)
      BB = 14.7D0*T1/(28.0D0+CAZ)
      AK = (LOG(AK)+CAZ*COS(AA)/(1.0D0+0.008D0*CAZ))/COS(BB)
      FK = 0.12125D0*AK*AK/CAZ + 1.5D0
  180 CONTINUE
C-----------------------------------------------------------------------
C     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
C-----------------------------------------------------------------------
      K = FK
      FK = K
      FKS = FK*FK
      P1R = CZEROR
      P1I = CZEROI
      P2R = TOL
      P2I = CZEROI
      CSR = P2R
      CSI = P2I
      DO 190 I=1,K
        A1 = FKS - FK
        AK = (FKS+FK)/(A1+FHS)
        RAK = 2.0D0/(FK+CONER)
        CBR = (FK+ZR)*RAK
        CBI = ZI*RAK
        PTR = P2R
        PTI = P2I
        P2R = (PTR*CBR-PTI*CBI-P1R)*AK
        P2I = (PTI*CBR+PTR*CBI-P1I)*AK
        P1R = PTR
        P1I = PTI
        CSR = CSR + P2R
        CSI = CSI + P2I
        FKS = A1 - FK + CONER
        FK = FK - CONER
  190 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE (P2/CS)=(P2/ABS(CS))*(CONJG(CS)/ABS(CS)) FOR BETTER
C     SCALING
C-----------------------------------------------------------------------
      TM = ZABS(CSR,CSI)
      PTR = 1.0D0/TM
      S1R = P2R*PTR
      S1I = P2I*PTR
      CSR = CSR*PTR
      CSI = -CSI*PTR
      CALL ZMLT(COEFR, COEFI, S1R, S1I, STR, STI)
      CALL ZMLT(STR, STI, CSR, CSI, S1R, S1I)
      IF (INU.GT.0 .OR. N.GT.1) GO TO 200
      ZDR = ZR
      ZDI = ZI
      IF(IFLAG.EQ.1) GO TO 270
      GO TO 240
  200 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE P1/P2=(P1/ABS(P2)*CONJG(P2)/ABS(P2) FOR SCALING
C-----------------------------------------------------------------------
      TM = ZABS(P2R,P2I)
      PTR = 1.0D0/TM
      P1R = P1R*PTR
      P1I = P1I*PTR
      P2R = P2R*PTR
      P2I = -P2I*PTR
      CALL ZMLT(P1R, P1I, P2R, P2I, PTR, PTI)
      STR = DNU + 0.5D0 - PTR
      STI = -PTI
      CALL ZDIV(STR, STI, ZR, ZI, STR, STI)
      STR = STR + 1.0D0
      CALL ZMLT(STR, STI, S1R, S1I, S2R, S2I)
C-----------------------------------------------------------------------
C     FORWARD RECURSION ON THE THREE TERM RECURSION WITH RELATION WITH
C     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
C-----------------------------------------------------------------------
  210 CONTINUE
      STR = DNU + 1.0D0
      CKR = STR*RZR
      CKI = STR*RZI
      IF (N.EQ.1) INU = INU - 1
      IF (INU.GT.0) GO TO 220
      IF (N.GT.1) GO TO 215
      S1R = S2R
      S1I = S2I
  215 CONTINUE
      ZDR = ZR
      ZDI = ZI
      IF(IFLAG.EQ.1) GO TO 270
      GO TO 240
  220 CONTINUE
      INUB = 1
      IF(IFLAG.EQ.1) GO TO 261
  225 CONTINUE
      P1R = CSRR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 230 I=INUB,INU
        STR = S2R
        STI = S2I
        S2R = CKR*STR - CKI*STI + S1R
        S2I = CKR*STI + CKI*STR + S1I
        S1R = STR
        S1I = STI
        CKR = CKR + RZR
        CKI = CKI + RZI
        IF (KFLAG.GE.3) GO TO 230
        P2R = S2R*P1R
        P2I = S2I*P1R
        STR = ABS(P2R)
        STI = ABS(P2I)
        P2M = MAX(STR,STI)
        IF (P2M.LE.ASCLE) GO TO 230
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1R = S1R*P1R
        S1I = S1I*P1R
        S2R = P2R
        S2I = P2I
        STR = CSSR(KFLAG)
        S1R = S1R*STR
        S1I = S1I*STR
        S2R = S2R*STR
        S2I = S2I*STR
        P1R = CSRR(KFLAG)
  230 CONTINUE
      IF (N.NE.1) GO TO 240
      S1R = S2R
      S1I = S2I
  240 CONTINUE
      STR = CSRR(KFLAG)
      YR(1) = S1R*STR
      YI(1) = S1I*STR
      IF (N.EQ.1) RETURN
      YR(2) = S2R*STR
      YI(2) = S2I*STR
      IF (N.EQ.2) RETURN
      KK = 2
  250 CONTINUE
      KK = KK + 1
      IF (KK.GT.N) RETURN
      P1R = CSRR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 260 I=KK,N
        P2R = S2R
        P2I = S2I
        S2R = CKR*P2R - CKI*P2I + S1R
        S2I = CKI*P2R + CKR*P2I + S1I
        S1R = P2R
        S1I = P2I
        CKR = CKR + RZR
        CKI = CKI + RZI
        P2R = S2R*P1R
        P2I = S2I*P1R
        YR(I) = P2R
        YI(I) = P2I
        IF (KFLAG.GE.3) GO TO 260
        STR = ABS(P2R)
        STI = ABS(P2I)
        P2M = MAX(STR,STI)
        IF (P2M.LE.ASCLE) GO TO 260
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1R = S1R*P1R
        S1I = S1I*P1R
        S2R = P2R
        S2I = P2I
        STR = CSSR(KFLAG)
        S1R = S1R*STR
        S1I = S1I*STR
        S2R = S2R*STR
        S2I = S2I*STR
        P1R = CSRR(KFLAG)
  260 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
C-----------------------------------------------------------------------
  261 CONTINUE
      HELIM = 0.5D0*ELIM
      ELM = EXP(-ELIM)
      CELMR = ELM
      ASCLE = BRY(1)
      ZDR = ZR
      ZDI = ZI
      IC = -1
      J = 2
      DO 262 I=1,INU
        STR = S2R
        STI = S2I
        S2R = STR*CKR-STI*CKI+S1R
        S2I = STI*CKR+STR*CKI+S1I
        S1R = STR
        S1I = STI
        CKR = CKR+RZR
        CKI = CKI+RZI
        AS = ZABS(S2R,S2I)
        ALAS = LOG(AS)
        P2R = -ZDR+ALAS
        IF(P2R.LT.(-ELIM)) GO TO 263
        CALL ZLOG(S2R,S2I,STR,STI,IDUM)
        P2R = -ZDR+STR
        P2I = -ZDI+STI
        P2M = EXP(P2R)/TOL
        P1R = P2M*COS(P2I)
        P1I = P2M*SIN(P2I)
        CALL ZUCHK(P1R,P1I,NW,ASCLE,TOL)
        IF(NW.NE.0) GO TO 263
        J = 3 - J
        CYR(J) = P1R
        CYI(J) = P1I
        IF(IC.EQ.(I-1)) GO TO 264
        IC = I
        GO TO 262
  263   CONTINUE
        IF(ALAS.LT.HELIM) GO TO 262
        ZDR = ZDR-ELIM
        S1R = S1R*CELMR
        S1I = S1I*CELMR
        S2R = S2R*CELMR
        S2I = S2I*CELMR
  262 CONTINUE
      IF(N.NE.1) GO TO 270
      S1R = S2R
      S1I = S2I
      GO TO 270
  264 CONTINUE
      KFLAG = 1
      INUB = I+1
      S2R = CYR(J)
      S2I = CYI(J)
      J = 3 - J
      S1R = CYR(J)
      S1I = CYI(J)
      IF(INUB.LE.INU) GO TO 225
      IF(N.NE.1) GO TO 240
      S1R = S2R
      S1I = S2I
      GO TO 240
  270 CONTINUE
      YR(1) = S1R
      YI(1) = S1I
      IF(N.EQ.1) GO TO 280
      YR(2) = S2R
      YI(2) = S2I
  280 CONTINUE
      ASCLE = BRY(1)
      CALL ZKSCL(ZDR,ZDI,FNU,N,YR,YI,NZ,RZR,RZI,ASCLE,TOL,ELIM)
      INU = N - NZ
      IF (INU.LE.0) RETURN
      KK = NZ + 1
      S1R = YR(KK)
      S1I = YI(KK)
      YR(KK) = S1R*CSRR(1)
      YI(KK) = S1I*CSRR(1)
      IF (INU.EQ.1) RETURN
      KK = NZ + 2
      S2R = YR(KK)
      S2I = YI(KK)
      YR(KK) = S2R*CSRR(1)
      YI(KK) = S2I*CSRR(1)
      IF (INU.EQ.2) RETURN
      T2 = FNU + (KK-1)
      CKR = T2*RZR
      CKI = T2*RZI
      KFLAG = 1
      GO TO 250
  290 CONTINUE
C-----------------------------------------------------------------------
C     SCALE BY EXP(Z), IFLAG = 1 CASES
C-----------------------------------------------------------------------
      KODED = 2
      IFLAG = 1
      KFLAG = 2
      GO TO 120
C-----------------------------------------------------------------------
C     FNU=HALF ODD INTEGER CASE, DNU=-0.5
C-----------------------------------------------------------------------
  300 CONTINUE
      S1R = COEFR
      S1I = COEFI
      S2R = COEFR
      S2I = COEFI
      GO TO 210
C
C
  310 CONTINUE
      NZ=-2
      RETURN
      END
*DECK ZBUNI
      SUBROUTINE ZBUNI (ZR, ZI, FNU, KODE, N, YR, YI, NZ, NUI, NLAST,
     +   FNUL, TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  ZBUNI
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CBUNI-A, ZBUNI-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZBUNI COMPUTES THE I BESSEL FUNCTION FOR LARGE ABS(Z).GT.
C     FNUL AND FNU+N-1.LT.FNUL. THE ORDER IS INCREASED FROM
C     FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING
C     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z)
C     ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2
C
C***SEE ALSO  ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZUNI1, ZUNI2
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZBUNI
C     COMPLEX CSCL,CSCR,CY,RZ,ST,S1,S2,Y,Z
      DOUBLE PRECISION ALIM, AX, AY, CSCLR, CSCRR, CYI, CYR, DFNU,
     * ELIM, FNU, FNUI, FNUL, GNU, RAZ, RZI, RZR, STI, STR, S1I, S1R,
     * S2I, S2R, TOL, YI, YR, ZI, ZR, ZABS, ASCLE, BRY, C1R, C1I, C1M,
     * D1MACH
      INTEGER I, IFLAG, IFORM, K, KODE, N, NL, NLAST, NUI, NW, NZ
      DIMENSION YR(N), YI(N), CYR(2), CYI(2), BRY(3)
      EXTERNAL ZABS
C***FIRST EXECUTABLE STATEMENT  ZBUNI
      NZ = 0
      AX = ABS(ZR)*1.7321D0
      AY = ABS(ZI)
      IFORM = 1
      IF (AY.GT.AX) IFORM = 2
      IF (NUI.EQ.0) GO TO 60
      FNUI = NUI
      DFNU = FNU + (N-1)
      GNU = DFNU + FNUI
      IF (IFORM.EQ.2) GO TO 10
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
      CALL ZUNI1(ZR, ZI, GNU, KODE, 2, CYR, CYI, NW, NLAST, FNUL, TOL,
     * ELIM, ALIM)
      GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
      CALL ZUNI2(ZR, ZI, GNU, KODE, 2, CYR, CYI, NW, NLAST, FNUL, TOL,
     * ELIM, ALIM)
   20 CONTINUE
      IF (NW.LT.0) GO TO 50
      IF (NW.NE.0) GO TO 90
      STR = ZABS(CYR(1),CYI(1))
C----------------------------------------------------------------------
C     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
C----------------------------------------------------------------------
      BRY(1)=1.0D+3*D1MACH(1)/TOL
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = BRY(2)
      IFLAG = 2
      ASCLE = BRY(2)
      CSCLR = 1.0D0
      IF (STR.GT.BRY(1)) GO TO 21
      IFLAG = 1
      ASCLE = BRY(1)
      CSCLR = 1.0D0/TOL
      GO TO 25
   21 CONTINUE
      IF (STR.LT.BRY(2)) GO TO 25
      IFLAG = 3
      ASCLE=BRY(3)
      CSCLR = TOL
   25 CONTINUE
      CSCRR = 1.0D0/CSCLR
      S1R = CYR(2)*CSCLR
      S1I = CYI(2)*CSCLR
      S2R = CYR(1)*CSCLR
      S2I = CYI(1)*CSCLR
      RAZ = 1.0D0/ZABS(ZR,ZI)
      STR = ZR*RAZ
      STI = -ZI*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      DO 30 I=1,NUI
        STR = S2R
        STI = S2I
        S2R = (DFNU+FNUI)*(RZR*STR-RZI*STI) + S1R
        S2I = (DFNU+FNUI)*(RZR*STI+RZI*STR) + S1I
        S1R = STR
        S1I = STI
        FNUI = FNUI - 1.0D0
        IF (IFLAG.GE.3) GO TO 30
        STR = S2R*CSCRR
        STI = S2I*CSCRR
        C1R = ABS(STR)
        C1I = ABS(STI)
        C1M = MAX(C1R,C1I)
        IF (C1M.LE.ASCLE) GO TO 30
        IFLAG = IFLAG+1
        ASCLE = BRY(IFLAG)
        S1R = S1R*CSCRR
        S1I = S1I*CSCRR
        S2R = STR
        S2I = STI
        CSCLR = CSCLR*TOL
        CSCRR = 1.0D0/CSCLR
        S1R = S1R*CSCLR
        S1I = S1I*CSCLR
        S2R = S2R*CSCLR
        S2I = S2I*CSCLR
   30 CONTINUE
      YR(N) = S2R*CSCRR
      YI(N) = S2I*CSCRR
      IF (N.EQ.1) RETURN
      NL = N - 1
      FNUI = NL
      K = NL
      DO 40 I=1,NL
        STR = S2R
        STI = S2I
        S2R = (FNU+FNUI)*(RZR*STR-RZI*STI) + S1R
        S2I = (FNU+FNUI)*(RZR*STI+RZI*STR) + S1I
        S1R = STR
        S1I = STI
        STR = S2R*CSCRR
        STI = S2I*CSCRR
        YR(K) = STR
        YI(K) = STI
        FNUI = FNUI - 1.0D0
        K = K - 1
        IF (IFLAG.GE.3) GO TO 40
        C1R = ABS(STR)
        C1I = ABS(STI)
        C1M = MAX(C1R,C1I)
        IF (C1M.LE.ASCLE) GO TO 40
        IFLAG = IFLAG+1
        ASCLE = BRY(IFLAG)
        S1R = S1R*CSCRR
        S1I = S1I*CSCRR
        S2R = STR
        S2I = STI
        CSCLR = CSCLR*TOL
        CSCRR = 1.0D0/CSCLR
        S1R = S1R*CSCLR
        S1I = S1I*CSCLR
        S2R = S2R*CSCLR
        S2I = S2I*CSCLR
   40 CONTINUE
      RETURN
   50 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
   60 CONTINUE
      IF (IFORM.EQ.2) GO TO 70
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
      CALL ZUNI1(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL,
     * ELIM, ALIM)
      GO TO 80
   70 CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
      CALL ZUNI2(ZR, ZI, FNU, KODE, N, YR, YI, NW, NLAST, FNUL, TOL,
     * ELIM, ALIM)
   80 CONTINUE
      IF (NW.LT.0) GO TO 50
      NZ = NW
      RETURN
   90 CONTINUE
      NLAST = N
      RETURN
      END
*DECK ZBUNK
      SUBROUTINE ZBUNK (ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
     +   ALIM)
C***BEGIN PROLOGUE  ZBUNK
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CBUNI-A, ZBUNI-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU.GT.FNUL.
C     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
C     IN ZUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN ZUNK2
C
C***SEE ALSO  ZBESH, ZBESK
C***ROUTINES CALLED  ZUNK1, ZUNK2
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZBUNK
C     COMPLEX Y,Z
      DOUBLE PRECISION ALIM, AX, AY, ELIM, FNU, TOL, YI, YR, ZI, ZR
      INTEGER KODE, MR, N, NZ
      DIMENSION YR(N), YI(N)
C***FIRST EXECUTABLE STATEMENT  ZBUNK
      NZ = 0
      AX = ABS(ZR)*1.7321D0
      AY = ABS(ZI)
      IF (AY.GT.AX) GO TO 10
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
C     -PI/3.LE.ARG(Z).LE.PI/3
C-----------------------------------------------------------------------
      CALL ZUNK1(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, ALIM)
      GO TO 20
   10 CONTINUE
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C     AND HPI=PI/2
C-----------------------------------------------------------------------
      CALL ZUNK2(ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM, ALIM)
   20 CONTINUE
      RETURN
      END
*DECK ZDIV
      SUBROUTINE ZDIV (AR, AI, BR, BI, CR, CI)
C***BEGIN PROLOGUE  ZDIV
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
C            ZBIRY
C***LIBRARY   SLATEC
C***TYPE      ALL (ZDIV-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     DOUBLE PRECISION COMPLEX DIVIDE C=A/B.
C
C***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
C***ROUTINES CALLED  ZABS
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZDIV
      DOUBLE PRECISION AR, AI, BR, BI, CR, CI, BM, CA, CB, CC, CD
      DOUBLE PRECISION ZABS
      EXTERNAL ZABS
C***FIRST EXECUTABLE STATEMENT  ZDIV
      BM = 1.0D0/ZABS(BR,BI)
      CC = BR*BM
      CD = BI*BM
      CA = (AR*CC+AI*CD)*BM
      CB = (AI*CC-AR*CD)*BM
      CR = CA
      CI = CB
      RETURN
      END
*DECK ZEXP
      SUBROUTINE ZEXP (AR, AI, BR, BI)
C***BEGIN PROLOGUE  ZEXP
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
C            ZBIRY
C***LIBRARY   SLATEC
C***TYPE      ALL (ZEXP-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     DOUBLE PRECISION COMPLEX EXPONENTIAL FUNCTION B=EXP(A)
C
C***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZEXP
      DOUBLE PRECISION AR, AI, BR, BI, ZM, CA, CB
C***FIRST EXECUTABLE STATEMENT  ZEXP
      ZM = EXP(AR)
      CA = ZM*COS(AI)
      CB = ZM*SIN(AI)
      BR = CA
      BI = CB
      RETURN
      END
*DECK ZKSCL
      SUBROUTINE ZKSCL (ZRR, ZRI, FNU, N, YR, YI, NZ, RZR, RZI, ASCLE,
     +   TOL, ELIM)
C***BEGIN PROLOGUE  ZKSCL
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CKSCL-A, ZKSCL-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE
C     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN
C     RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.
C
C***SEE ALSO  ZBESK
C***ROUTINES CALLED  ZABS, ZLOG, ZUCHK
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   930122  Added ZLOG to EXTERNAL statement.  (RWC)
C***END PROLOGUE  ZKSCL
C     COMPLEX CK,CS,CY,CZERO,RZ,S1,S2,Y,ZR,ZD,CELM
      DOUBLE PRECISION ACS, AS, ASCLE, CKI, CKR, CSI, CSR, CYI,
     * CYR, ELIM, FN, FNU, RZI, RZR, STR, S1I, S1R, S2I,
     * S2R, TOL, YI, YR, ZEROI, ZEROR, ZRI, ZRR, ZABS,
     * ZDR, ZDI, CELMR, ELM, HELIM, ALAS
      INTEGER I, IC, IDUM, KK, N, NN, NW, NZ
      DIMENSION YR(N), YI(N), CYR(2), CYI(2)
      EXTERNAL ZABS, ZLOG
      DATA ZEROR,ZEROI / 0.0D0 , 0.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZKSCL
      NZ = 0
      IC = 0
      NN = MIN(2,N)
      DO 10 I=1,NN
        S1R = YR(I)
        S1I = YI(I)
        CYR(I) = S1R
        CYI(I) = S1I
        AS = ZABS(S1R,S1I)
        ACS = -ZRR + LOG(AS)
        NZ = NZ + 1
        YR(I) = ZEROR
        YI(I) = ZEROI
        IF (ACS.LT.(-ELIM)) GO TO 10
        CALL ZLOG(S1R, S1I, CSR, CSI, IDUM)
        CSR = CSR - ZRR
        CSI = CSI - ZRI
        STR = EXP(CSR)/TOL
        CSR = STR*COS(CSI)
        CSI = STR*SIN(CSI)
        CALL ZUCHK(CSR, CSI, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 10
        YR(I) = CSR
        YI(I) = CSI
        IC = I
        NZ = NZ - 1
   10 CONTINUE
      IF (N.EQ.1) RETURN
      IF (IC.GT.1) GO TO 20
      YR(1) = ZEROR
      YI(1) = ZEROI
      NZ = 2
   20 CONTINUE
      IF (N.EQ.2) RETURN
      IF (NZ.EQ.0) RETURN
      FN = FNU + 1.0D0
      CKR = FN*RZR
      CKI = FN*RZI
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      HELIM = 0.5D0*ELIM
      ELM = EXP(-ELIM)
      CELMR = ELM
      ZDR = ZRR
      ZDI = ZRI
C
C     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF
C     S2 GETS LARGER THAN EXP(ELIM/2)
C
      DO 30 I=3,N
        KK = I
        CSR = S2R
        CSI = S2I
        S2R = CKR*CSR - CKI*CSI + S1R
        S2I = CKI*CSR + CKR*CSI + S1I
        S1R = CSR
        S1I = CSI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AS = ZABS(S2R,S2I)
        ALAS = LOG(AS)
        ACS = -ZDR + ALAS
        NZ = NZ + 1
        YR(I) = ZEROR
        YI(I) = ZEROI
        IF (ACS.LT.(-ELIM)) GO TO 25
        CALL ZLOG(S2R, S2I, CSR, CSI, IDUM)
        CSR = CSR - ZDR
        CSI = CSI - ZDI
        STR = EXP(CSR)/TOL
        CSR = STR*COS(CSI)
        CSI = STR*SIN(CSI)
        CALL ZUCHK(CSR, CSI, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 25
        YR(I) = CSR
        YI(I) = CSI
        NZ = NZ - 1
        IF (IC.EQ.KK-1) GO TO 40
        IC = KK
        GO TO 30
   25   CONTINUE
        IF(ALAS.LT.HELIM) GO TO 30
        ZDR = ZDR - ELIM
        S1R = S1R*CELMR
        S1I = S1I*CELMR
        S2R = S2R*CELMR
        S2I = S2I*CELMR
   30 CONTINUE
      NZ = N
      IF(IC.EQ.N) NZ=N-1
      GO TO 45
   40 CONTINUE
      NZ = KK - 2
   45 CONTINUE
      DO 50 I=1,NZ
        YR(I) = ZEROR
        YI(I) = ZEROI
   50 CONTINUE
      RETURN
      END
*DECK ZLOG
      SUBROUTINE ZLOG (AR, AI, BR, BI, IERR)
C***BEGIN PROLOGUE  ZLOG
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
C            ZBIRY
C***LIBRARY   SLATEC
C***TYPE      ALL (ZLOG-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     DOUBLE PRECISION COMPLEX LOGARITHM B=CLOG(A)
C     IERR=0,NORMAL RETURN      IERR=1, Z=CMPLX(0.0,0.0)
C***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
C***ROUTINES CALLED  ZABS
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZLOG
      DOUBLE PRECISION AR, AI, BR, BI, ZM, DTHETA, DPI, DHPI
      DOUBLE PRECISION ZABS
      INTEGER IERR
      EXTERNAL ZABS
      DATA DPI , DHPI  / 3.141592653589793238462643383D+0,
     1                   1.570796326794896619231321696D+0/
C***FIRST EXECUTABLE STATEMENT  ZLOG
      IERR=0
      IF (AR.EQ.0.0D+0) GO TO 10
      IF (AI.EQ.0.0D+0) GO TO 20
      DTHETA = DATAN(AI/AR)
      IF (DTHETA.LE.0.0D+0) GO TO 40
      IF (AR.LT.0.0D+0) DTHETA = DTHETA - DPI
      GO TO 50
   10 IF (AI.EQ.0.0D+0) GO TO 60
      BI = DHPI
      BR = LOG(ABS(AI))
      IF (AI.LT.0.0D+0) BI = -BI
      RETURN
   20 IF (AR.GT.0.0D+0) GO TO 30
      BR = LOG(ABS(AR))
      BI = DPI
      RETURN
   30 BR = LOG(AR)
      BI = 0.0D+0
      RETURN
   40 IF (AR.LT.0.0D+0) DTHETA = DTHETA + DPI
   50 ZM = ZABS(AR,AI)
      BR = LOG(ZM)
      BI = DTHETA
      RETURN
   60 CONTINUE
      IERR=1
      RETURN
      END
*DECK ZMLRI
      SUBROUTINE ZMLRI (ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL)
C***BEGIN PROLOGUE  ZMLRI
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CMLRI-A, ZMLRI-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY THE
C     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
C
C***SEE ALSO  ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, DGAMLN, ZABS, ZEXP, ZLOG, ZMLT
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   930122  Added ZEXP and ZLOG to EXTERNAL statement.  (RWC)
C***END PROLOGUE  ZMLRI
C     COMPLEX CK,CNORM,CONE,CTWO,CZERO,PT,P1,P2,RZ,SUM,Y,Z
      DOUBLE PRECISION ACK, AK, AP, AT, AZ, BK, CKI, CKR, CNORMI,
     * CNORMR, CONEI, CONER, FKAP, FKK, FLAM, FNF, FNU, PTI, PTR, P1I,
     * P1R, P2I, P2R, RAZ, RHO, RHO2, RZI, RZR, SCLE, STI, STR, SUMI,
     * SUMR, TFNF, TOL, TST, YI, YR, ZEROI, ZEROR, ZI, ZR, DGAMLN,
     * D1MACH, ZABS
      INTEGER I, IAZ, IDUM, IFNU, INU, ITIME, K, KK, KM, KODE, M, N, NZ
      DIMENSION YR(N), YI(N)
      EXTERNAL ZABS, ZEXP, ZLOG
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZMLRI
      SCLE = D1MACH(1)/TOL
      NZ=0
      AZ = ZABS(ZR,ZI)
      IAZ = AZ
      IFNU = FNU
      INU = IFNU + N - 1
      AT = IAZ + 1.0D0
      RAZ = 1.0D0/AZ
      STR = ZR*RAZ
      STI = -ZI*RAZ
      CKR = STR*AT*RAZ
      CKI = STI*AT*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      P1R = ZEROR
      P1I = ZEROI
      P2R = CONER
      P2I = CONEI
      ACK = (AT+1.0D0)*RAZ
      RHO = ACK + SQRT(ACK*ACK-1.0D0)
      RHO2 = RHO*RHO
      TST = (RHO2+RHO2)/((RHO2-1.0D0)*(RHO-1.0D0))
      TST = TST/TOL
C-----------------------------------------------------------------------
C     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
C-----------------------------------------------------------------------
      AK = AT
      DO 10 I=1,80
        PTR = P2R
        PTI = P2I
        P2R = P1R - (CKR*PTR-CKI*PTI)
        P2I = P1I - (CKI*PTR+CKR*PTI)
        P1R = PTR
        P1I = PTI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AP = ZABS(P2R,P2I)
        IF (AP.GT.TST*AK*AK) GO TO 20
        AK = AK + 1.0D0
   10 CONTINUE
      GO TO 110
   20 CONTINUE
      I = I + 1
      K = 0
      IF (INU.LT.IAZ) GO TO 40
C-----------------------------------------------------------------------
C     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
C-----------------------------------------------------------------------
      P1R = ZEROR
      P1I = ZEROI
      P2R = CONER
      P2I = CONEI
      AT = INU + 1.0D0
      STR = ZR*RAZ
      STI = -ZI*RAZ
      CKR = STR*AT*RAZ
      CKI = STI*AT*RAZ
      ACK = AT*RAZ
      TST = SQRT(ACK/TOL)
      ITIME = 1
      DO 30 K=1,80
        PTR = P2R
        PTI = P2I
        P2R = P1R - (CKR*PTR-CKI*PTI)
        P2I = P1I - (CKR*PTI+CKI*PTR)
        P1R = PTR
        P1I = PTI
        CKR = CKR + RZR
        CKI = CKI + RZI
        AP = ZABS(P2R,P2I)
        IF (AP.LT.TST) GO TO 30
        IF (ITIME.EQ.2) GO TO 40
        ACK = ZABS(CKR,CKI)
        FLAM = ACK + SQRT(ACK*ACK-1.0D0)
        FKAP = AP/ZABS(P1R,P1I)
        RHO = MIN(FLAM,FKAP)
        TST = TST*SQRT(RHO/(RHO*RHO-1.0D0))
        ITIME = 2
   30 CONTINUE
      GO TO 110
   40 CONTINUE
C-----------------------------------------------------------------------
C     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
C-----------------------------------------------------------------------
      K = K + 1
      KK = MAX(I+IAZ,K+INU)
      FKK = KK
      P1R = ZEROR
      P1I = ZEROI
C-----------------------------------------------------------------------
C     SCALE P2 AND SUM BY SCLE
C-----------------------------------------------------------------------
      P2R = SCLE
      P2I = ZEROI
      FNF = FNU - IFNU
      TFNF = FNF + FNF
      BK = DGAMLN(FKK+TFNF+1.0D0,IDUM) - DGAMLN(FKK+1.0D0,IDUM) -
     * DGAMLN(TFNF+1.0D0,IDUM)
      BK = EXP(BK)
      SUMR = ZEROR
      SUMI = ZEROI
      KM = KK - INU
      DO 50 I=1,KM
        PTR = P2R
        PTI = P2I
        P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
        P2I = P1I + (FKK+FNF)*(RZI*PTR+RZR*PTI)
        P1R = PTR
        P1I = PTI
        AK = 1.0D0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUMR = SUMR + (ACK+BK)*P1R
        SUMI = SUMI + (ACK+BK)*P1I
        BK = ACK
        FKK = FKK - 1.0D0
   50 CONTINUE
      YR(N) = P2R
      YI(N) = P2I
      IF (N.EQ.1) GO TO 70
      DO 60 I=2,N
        PTR = P2R
        PTI = P2I
        P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
        P2I = P1I + (FKK+FNF)*(RZI*PTR+RZR*PTI)
        P1R = PTR
        P1I = PTI
        AK = 1.0D0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUMR = SUMR + (ACK+BK)*P1R
        SUMI = SUMI + (ACK+BK)*P1I
        BK = ACK
        FKK = FKK - 1.0D0
        M = N - I + 1
        YR(M) = P2R
        YI(M) = P2I
   60 CONTINUE
   70 CONTINUE
      IF (IFNU.LE.0) GO TO 90
      DO 80 I=1,IFNU
        PTR = P2R
        PTI = P2I
        P2R = P1R + (FKK+FNF)*(RZR*PTR-RZI*PTI)
        P2I = P1I + (FKK+FNF)*(RZR*PTI+RZI*PTR)
        P1R = PTR
        P1I = PTI
        AK = 1.0D0 - TFNF/(FKK+TFNF)
        ACK = BK*AK
        SUMR = SUMR + (ACK+BK)*P1R
        SUMI = SUMI + (ACK+BK)*P1I
        BK = ACK
        FKK = FKK - 1.0D0
   80 CONTINUE
   90 CONTINUE
      PTR = ZR
      PTI = ZI
      IF (KODE.EQ.2) PTR = ZEROR
      CALL ZLOG(RZR, RZI, STR, STI, IDUM)
      P1R = -FNF*STR + PTR
      P1I = -FNF*STI + PTI
      AP = DGAMLN(1.0D0+FNF,IDUM)
      PTR = P1R - AP
      PTI = P1I
C-----------------------------------------------------------------------
C     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
C     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
C-----------------------------------------------------------------------
      P2R = P2R + SUMR
      P2I = P2I + SUMI
      AP = ZABS(P2R,P2I)
      P1R = 1.0D0/AP
      CALL ZEXP(PTR, PTI, STR, STI)
      CKR = STR*P1R
      CKI = STI*P1R
      PTR = P2R*P1R
      PTI = -P2I*P1R
      CALL ZMLT(CKR, CKI, PTR, PTI, CNORMR, CNORMI)
      DO 100 I=1,N
        STR = YR(I)*CNORMR - YI(I)*CNORMI
        YI(I) = YR(I)*CNORMI + YI(I)*CNORMR
        YR(I) = STR
  100 CONTINUE
      RETURN
  110 CONTINUE
      NZ=-2
      RETURN
      END
*DECK ZMLT
      SUBROUTINE ZMLT (AR, AI, BR, BI, CR, CI)
C***BEGIN PROLOGUE  ZMLT
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
C            ZBIRY
C***LIBRARY   SLATEC
C***TYPE      ALL (ZMLT-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     DOUBLE PRECISION COMPLEX MULTIPLY, C=A*B.
C
C***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZMLT
      DOUBLE PRECISION AR, AI, BR, BI, CR, CI, CA, CB
C***FIRST EXECUTABLE STATEMENT  ZMLT
      CA = AR*BR - AI*BI
      CB = AR*BI + AI*BR
      CR = CA
      CI = CB
      RETURN
      END
*DECK ZRATI
      SUBROUTINE ZRATI (ZR, ZI, FNU, N, CYR, CYI, TOL)
C***BEGIN PROLOGUE  ZRATI
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH, ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CRATI-A, ZRATI-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZRATI COMPUTES RATIOS OF I BESSEL FUNCTIONS BY BACKWARD
C     RECURRENCE.  THE STARTING INDEX IS DETERMINED BY FORWARD
C     RECURRENCE AS DESCRIBED IN J. RES. OF NAT. BUR. OF STANDARDS-B,
C     MATHEMATICAL SCIENCES, VOL 77B, P111-114, SEPTEMBER, 1973,
C     BESSEL FUNCTIONS I AND J OF COMPLEX ARGUMENT AND INTEGER ORDER,
C     BY D. J. SOOKNE.
C
C***SEE ALSO  ZBESH, ZBESI, ZBESK
C***ROUTINES CALLED  ZABS, ZDIV
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZRATI
      DOUBLE PRECISION AK, AMAGZ, AP1, AP2, ARG, AZ, CDFNUI, CDFNUR,
     * CONEI, CONER, CYI, CYR, CZEROI, CZEROR, DFNU, FDNU, FLAM, FNU,
     * FNUP, PTI, PTR, P1I, P1R, P2I, P2R, RAK, RAP1, RHO, RT2, RZI,
     * RZR, TEST, TEST1, TOL, TTI, TTR, T1I, T1R, ZI, ZR, ZABS
      INTEGER I, ID, IDNU, INU, ITIME, K, KK, MAGZ, N
      DIMENSION CYR(N), CYI(N)
      EXTERNAL ZABS
      DATA CZEROR,CZEROI,CONER,CONEI,RT2/
     1 0.0D0, 0.0D0, 1.0D0, 0.0D0, 1.41421356237309505D0 /
C***FIRST EXECUTABLE STATEMENT  ZRATI
      AZ = ZABS(ZR,ZI)
      INU = FNU
      IDNU = INU + N - 1
      MAGZ = AZ
      AMAGZ = MAGZ+1
      FDNU = IDNU
      FNUP = MAX(AMAGZ,FDNU)
      ID = IDNU - MAGZ - 1
      ITIME = 1
      K = 1
      PTR = 1.0D0/AZ
      RZR = PTR*(ZR+ZR)*PTR
      RZI = -PTR*(ZI+ZI)*PTR
      T1R = RZR*FNUP
      T1I = RZI*FNUP
      P2R = -T1R
      P2I = -T1I
      P1R = CONER
      P1I = CONEI
      T1R = T1R + RZR
      T1I = T1I + RZI
      IF (ID.GT.0) ID = 0
      AP2 = ZABS(P2R,P2I)
      AP1 = ZABS(P1R,P1I)
C-----------------------------------------------------------------------
C     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNU
C     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
C     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
C     PREMATURELY.
C-----------------------------------------------------------------------
      ARG = (AP2+AP2)/(AP1*TOL)
      TEST1 = SQRT(ARG)
      TEST = TEST1
      RAP1 = 1.0D0/AP1
      P1R = P1R*RAP1
      P1I = P1I*RAP1
      P2R = P2R*RAP1
      P2I = P2I*RAP1
      AP2 = AP2*RAP1
   10 CONTINUE
      K = K + 1
      AP1 = AP2
      PTR = P2R
      PTI = P2I
      P2R = P1R - (T1R*PTR-T1I*PTI)
      P2I = P1I - (T1R*PTI+T1I*PTR)
      P1R = PTR
      P1I = PTI
      T1R = T1R + RZR
      T1I = T1I + RZI
      AP2 = ZABS(P2R,P2I)
      IF (AP1.LE.TEST) GO TO 10
      IF (ITIME.EQ.2) GO TO 20
      AK = ZABS(T1R,T1I)*0.5D0
      FLAM = AK + SQRT(AK*AK-1.0D0)
      RHO = MIN(AP2/AP1,FLAM)
      TEST = TEST1*SQRT(RHO/(RHO*RHO-1.0D0))
      ITIME = 2
      GO TO 10
   20 CONTINUE
      KK = K + 1 - ID
      AK = KK
      T1R = AK
      T1I = CZEROI
      DFNU = FNU + (N-1)
      P1R = 1.0D0/AP2
      P1I = CZEROI
      P2R = CZEROR
      P2I = CZEROI
      DO 30 I=1,KK
        PTR = P1R
        PTI = P1I
        RAP1 = DFNU + T1R
        TTR = RZR*RAP1
        TTI = RZI*RAP1
        P1R = (PTR*TTR-PTI*TTI) + P2R
        P1I = (PTR*TTI+PTI*TTR) + P2I
        P2R = PTR
        P2I = PTI
        T1R = T1R - CONER
   30 CONTINUE
      IF (P1R.NE.CZEROR .OR. P1I.NE.CZEROI) GO TO 40
      P1R = TOL
      P1I = TOL
   40 CONTINUE
      CALL ZDIV(P2R, P2I, P1R, P1I, CYR(N), CYI(N))
      IF (N.EQ.1) RETURN
      K = N - 1
      AK = K
      T1R = AK
      T1I = CZEROI
      CDFNUR = FNU*RZR
      CDFNUI = FNU*RZI
      DO 60 I=2,N
        PTR = CDFNUR + (T1R*RZR-T1I*RZI) + CYR(K+1)
        PTI = CDFNUI + (T1R*RZI+T1I*RZR) + CYI(K+1)
        AK = ZABS(PTR,PTI)
        IF (AK.NE.CZEROR) GO TO 50
        PTR = TOL
        PTI = TOL
        AK = TOL*RT2
   50   CONTINUE
        RAK = CONER/AK
        CYR(K) = RAK*PTR*RAK
        CYI(K) = -RAK*PTI*RAK
        T1R = T1R - CONER
        K = K - 1
   60 CONTINUE
      RETURN
      END
*DECK ZS1S2
      SUBROUTINE ZS1S2 (ZRR, ZRI, S1R, S1I, S2R, S2I, NZ, ASCLE, ALIM,
     +   IUF)
C***BEGIN PROLOGUE  ZS1S2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZAIRY and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CS1S2-A, ZS1S2-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE
C     ADDITION OF THE I AND K FUNCTIONS IN THE ANALYTIC CON-
C     TINUATION FORMULA WHERE S1=K FUNCTION AND S2=I FUNCTION.
C     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF
C     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER
C     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE
C     PRECISION ABOVE THE UNDERFLOW LIMIT.
C
C***SEE ALSO  ZAIRY, ZBESK
C***ROUTINES CALLED  ZABS, ZEXP, ZLOG
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   930122  Added ZEXP and ZLOG to EXTERNAL statement.  (RWC)
C***END PROLOGUE  ZS1S2
C     COMPLEX CZERO,C1,S1,S1D,S2,ZR
      DOUBLE PRECISION AA, ALIM, ALN, ASCLE, AS1, AS2, C1I, C1R, S1DI,
     * S1DR, S1I, S1R, S2I, S2R, ZEROI, ZEROR, ZRI, ZRR, ZABS
      INTEGER IUF, IDUM, NZ
      EXTERNAL ZABS, ZEXP, ZLOG
      DATA ZEROR,ZEROI  / 0.0D0 , 0.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZS1S2
      NZ = 0
      AS1 = ZABS(S1R,S1I)
      AS2 = ZABS(S2R,S2I)
      IF (S1R.EQ.0.0D0 .AND. S1I.EQ.0.0D0) GO TO 10
      IF (AS1.EQ.0.0D0) GO TO 10
      ALN = -ZRR - ZRR + LOG(AS1)
      S1DR = S1R
      S1DI = S1I
      S1R = ZEROR
      S1I = ZEROI
      AS1 = ZEROR
      IF (ALN.LT.(-ALIM)) GO TO 10
      CALL ZLOG(S1DR, S1DI, C1R, C1I, IDUM)
      C1R = C1R - ZRR - ZRR
      C1I = C1I - ZRI - ZRI
      CALL ZEXP(C1R, C1I, S1R, S1I)
      AS1 = ZABS(S1R,S1I)
      IUF = IUF + 1
   10 CONTINUE
      AA = MAX(AS1,AS2)
      IF (AA.GT.ASCLE) RETURN
      S1R = ZEROR
      S1I = ZEROI
      S2R = ZEROR
      S2I = ZEROI
      NZ = 1
      IUF = 0
      RETURN
      END
*DECK ZSERI
      SUBROUTINE ZSERI (ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM,
     +   ALIM)
C***BEGIN PROLOGUE  ZSERI
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CSERI-A, ZSERI-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
C     MEANS OF THE POWER SERIES FOR LARGE ABS(Z) IN THE
C     REGION ABS(Z).LE.2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
C     NZ.GT.0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
C     DUE TO UNDERFLOW. NZ.LT.0 MEANS UNDERFLOW OCCURRED, BUT THE
C     CONDITION ABS(Z).LE.2*SQRT(FNU+1) WAS VIOLATED AND THE
C     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
C
C***SEE ALSO  ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, DGAMLN, ZABS, ZDIV, ZLOG, ZMLT, ZUCHK
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   930122  Added ZLOG to EXTERNAL statement.  (RWC)
C***END PROLOGUE  ZSERI
C     COMPLEX AK1,CK,COEF,CONE,CRSC,CSCL,CZ,CZERO,HZ,RZ,S1,S2,Y,Z
      DOUBLE PRECISION AA, ACZ, AK, AK1I, AK1R, ALIM, ARM, ASCLE, ATOL,
     * AZ, CKI, CKR, COEFI, COEFR, CONEI, CONER, CRSCR, CZI, CZR, DFNU,
     * ELIM, FNU, FNUP, HZI, HZR, RAZ, RS, RTR1, RZI, RZR, S, SS, STI,
     * STR, S1I, S1R, S2I, S2R, TOL, YI, YR, WI, WR, ZEROI, ZEROR, ZI,
     * ZR, DGAMLN, D1MACH, ZABS
      INTEGER I, IB, IDUM, IFLAG, IL, K, KODE, L, M, N, NN, NZ, NW
      DIMENSION YR(N), YI(N), WR(2), WI(2)
      EXTERNAL ZABS, ZLOG
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZSERI
      NZ = 0
      AZ = ZABS(ZR,ZI)
      IF (AZ.EQ.0.0D0) GO TO 160
      ARM = 1.0D+3*D1MACH(1)
      RTR1 = SQRT(ARM)
      CRSCR = 1.0D0
      IFLAG = 0
      IF (AZ.LT.ARM) GO TO 150
      HZR = 0.5D0*ZR
      HZI = 0.5D0*ZI
      CZR = ZEROR
      CZI = ZEROI
      IF (AZ.LE.RTR1) GO TO 10
      CALL ZMLT(HZR, HZI, HZR, HZI, CZR, CZI)
   10 CONTINUE
      ACZ = ZABS(CZR,CZI)
      NN = N
      CALL ZLOG(HZR, HZI, CKR, CKI, IDUM)
   20 CONTINUE
      DFNU = FNU + (NN-1)
      FNUP = DFNU + 1.0D0
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
      AK1R = CKR*DFNU
      AK1I = CKI*DFNU
      AK = DGAMLN(FNUP,IDUM)
      AK1R = AK1R - AK
      IF (KODE.EQ.2) AK1R = AK1R - ZR
      IF (AK1R.GT.(-ELIM)) GO TO 40
   30 CONTINUE
      NZ = NZ + 1
      YR(NN) = ZEROR
      YI(NN) = ZEROI
      IF (ACZ.GT.DFNU) GO TO 190
      NN = NN - 1
      IF (NN.EQ.0) RETURN
      GO TO 20
   40 CONTINUE
      IF (AK1R.GT.(-ALIM)) GO TO 50
      IFLAG = 1
      SS = 1.0D0/TOL
      CRSCR = TOL
      ASCLE = ARM*SS
   50 CONTINUE
      AA = EXP(AK1R)
      IF (IFLAG.EQ.1) AA = AA*SS
      COEFR = AA*COS(AK1I)
      COEFI = AA*SIN(AK1I)
      ATOL = TOL*ACZ/FNUP
      IL = MIN(2,NN)
      DO 90 I=1,IL
        DFNU = FNU + (NN-I)
        FNUP = DFNU + 1.0D0
        S1R = CONER
        S1I = CONEI
        IF (ACZ.LT.TOL*FNUP) GO TO 70
        AK1R = CONER
        AK1I = CONEI
        AK = FNUP + 2.0D0
        S = FNUP
        AA = 2.0D0
   60   CONTINUE
        RS = 1.0D0/S
        STR = AK1R*CZR - AK1I*CZI
        STI = AK1R*CZI + AK1I*CZR
        AK1R = STR*RS
        AK1I = STI*RS
        S1R = S1R + AK1R
        S1I = S1I + AK1I
        S = S + AK
        AK = AK + 2.0D0
        AA = AA*ACZ*RS
        IF (AA.GT.ATOL) GO TO 60
   70   CONTINUE
        S2R = S1R*COEFR - S1I*COEFI
        S2I = S1R*COEFI + S1I*COEFR
        WR(I) = S2R
        WI(I) = S2I
        IF (IFLAG.EQ.0) GO TO 80
        CALL ZUCHK(S2R, S2I, NW, ASCLE, TOL)
        IF (NW.NE.0) GO TO 30
   80   CONTINUE
        M = NN - I + 1
        YR(M) = S2R*CRSCR
        YI(M) = S2I*CRSCR
        IF (I.EQ.IL) GO TO 90
        CALL ZDIV(COEFR, COEFI, HZR, HZI, STR, STI)
        COEFR = STR*DFNU
        COEFI = STI*DFNU
   90 CONTINUE
      IF (NN.LE.2) RETURN
      K = NN - 2
      AK = K
      RAZ = 1.0D0/AZ
      STR = ZR*RAZ
      STI = -ZI*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      IF (IFLAG.EQ.1) GO TO 120
      IB = 3
  100 CONTINUE
      DO 110 I=IB,NN
        YR(K) = (AK+FNU)*(RZR*YR(K+1)-RZI*YI(K+1)) + YR(K+2)
        YI(K) = (AK+FNU)*(RZR*YI(K+1)+RZI*YR(K+1)) + YI(K+2)
        AK = AK - 1.0D0
        K = K - 1
  110 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD WITH SCALED VALUES
C-----------------------------------------------------------------------
  120 CONTINUE
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE
C     UNDERFLOW LIMIT = ASCLE = D1MACH(1)*SS*1.0D+3
C-----------------------------------------------------------------------
      S1R = WR(1)
      S1I = WI(1)
      S2R = WR(2)
      S2I = WI(2)
      DO 130 L=3,NN
        CKR = S2R
        CKI = S2I
        S2R = S1R + (AK+FNU)*(RZR*CKR-RZI*CKI)
        S2I = S1I + (AK+FNU)*(RZR*CKI+RZI*CKR)
        S1R = CKR
        S1I = CKI
        CKR = S2R*CRSCR
        CKI = S2I*CRSCR
        YR(K) = CKR
        YI(K) = CKI
        AK = AK - 1.0D0
        K = K - 1
        IF (ZABS(CKR,CKI).GT.ASCLE) GO TO 140
  130 CONTINUE
      RETURN
  140 CONTINUE
      IB = L + 1
      IF (IB.GT.NN) RETURN
      GO TO 100
  150 CONTINUE
      NZ = N
      IF (FNU.EQ.0.0D0) NZ = NZ - 1
  160 CONTINUE
      YR(1) = ZEROR
      YI(1) = ZEROI
      IF (FNU.NE.0.0D0) GO TO 170
      YR(1) = CONER
      YI(1) = CONEI
  170 CONTINUE
      IF (N.EQ.1) RETURN
      DO 180 I=2,N
        YR(I) = ZEROR
        YI(I) = ZEROI
  180 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     RETURN WITH NZ.LT.0 IF ABS(Z*Z/4).GT.FNU+N-NZ-1 COMPLETE
C     THE CALCULATION IN CBINU WITH N=N-ABS(NZ)
C-----------------------------------------------------------------------
  190 CONTINUE
      NZ = -NZ
      RETURN
      END
*DECK ZSHCH
      SUBROUTINE ZSHCH (ZR, ZI, CSHR, CSHI, CCHR, CCHI)
C***BEGIN PROLOGUE  ZSHCH
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CSHCH-A, ZSHCH-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y)
C     AND CCH=COSH(X+I*Y), WHERE I**2=-1.
C
C***SEE ALSO  ZBESH, ZBESK
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZSHCH
C
      DOUBLE PRECISION CCHI, CCHR, CH, CN, CSHI, CSHR, SH, SN, ZI, ZR
C***FIRST EXECUTABLE STATEMENT  ZSHCH
      SH = SINH(ZR)
      CH = COSH(ZR)
      SN = SIN(ZI)
      CN = COS(ZI)
      CSHR = SH*CN
      CSHI = CH*SN
      CCHR = CH*CN
      CCHI = SH*SN
      RETURN
      END
*DECK ZSQRT
      SUBROUTINE ZSQRT (AR, AI, BR, BI)
C***BEGIN PROLOGUE  ZSQRT
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
C            ZBIRY
C***LIBRARY   SLATEC
C***TYPE      ALL (ZSQRT-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     DOUBLE PRECISION COMPLEX SQUARE ROOT, B=CSQRT(A)
C
C***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
C***ROUTINES CALLED  ZABS
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZSQRT
      DOUBLE PRECISION AR, AI, BR, BI, ZM, DTHETA, DPI, DRT
      DOUBLE PRECISION ZABS
      EXTERNAL ZABS
      DATA DRT , DPI / 7.071067811865475244008443621D-1,
     1                 3.141592653589793238462643383D+0/
C***FIRST EXECUTABLE STATEMENT  ZSQRT
      ZM = ZABS(AR,AI)
      ZM = SQRT(ZM)
      IF (AR.EQ.0.0D+0) GO TO 10
      IF (AI.EQ.0.0D+0) GO TO 20
      DTHETA = DATAN(AI/AR)
      IF (DTHETA.LE.0.0D+0) GO TO 40
      IF (AR.LT.0.0D+0) DTHETA = DTHETA - DPI
      GO TO 50
   10 IF (AI.GT.0.0D+0) GO TO 60
      IF (AI.LT.0.0D+0) GO TO 70
      BR = 0.0D+0
      BI = 0.0D+0
      RETURN
   20 IF (AR.GT.0.0D+0) GO TO 30
      BR = 0.0D+0
      BI = SQRT(ABS(AR))
      RETURN
   30 BR = SQRT(AR)
      BI = 0.0D+0
      RETURN
   40 IF (AR.LT.0.0D+0) DTHETA = DTHETA + DPI
   50 DTHETA = DTHETA*0.5D+0
      BR = ZM*COS(DTHETA)
      BI = ZM*SIN(DTHETA)
      RETURN
   60 BR = ZM*DRT
      BI = ZM*DRT
      RETURN
   70 BR = ZM*DRT
      BI = -ZM*DRT
      RETURN
      END
*DECK ZUCHK
      SUBROUTINE ZUCHK (YR, YI, NZ, ASCLE, TOL)
C***BEGIN PROLOGUE  ZUCHK
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SERI, ZUOIK, ZUNK1, ZUNK2, ZUNI1, ZUNI2 and
C            ZKSCL
C***LIBRARY   SLATEC
C***TYPE      ALL (CUCHK-A, ZUCHK-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C      Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN
C      EXP(-ALIM)=ASCLE=1.0E+3*D1MACH(1)/TOL. THE TEST IS MADE TO SEE
C      IF THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDERFLOW
C      WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED
C      IF THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE
C      OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE
C      ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED.
C
C***SEE ALSO  SERI, ZKSCL, ZUNI1, ZUNI2, ZUNK1, ZUNK2, ZUOIK
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZUCHK
C
C     COMPLEX Y
      DOUBLE PRECISION ASCLE, SS, ST, TOL, WR, WI, YR, YI
      INTEGER NZ
C***FIRST EXECUTABLE STATEMENT  ZUCHK
      NZ = 0
      WR = ABS(YR)
      WI = ABS(YI)
      ST = MIN(WR,WI)
      IF (ST.GT.ASCLE) RETURN
      SS = MAX(WR,WI)
      ST = ST/TOL
      IF (SS.LT.ST) NZ = 1
      RETURN
      END
*DECK ZUNHJ
      SUBROUTINE ZUNHJ (ZR, ZI, FNU, IPMTR, TOL, PHIR, PHII, ARGR, ARGI,
     +   ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
C***BEGIN PROLOGUE  ZUNHJ
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CUNHJ-A, ZUNHJ-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     REFERENCES
C         HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND I.A.
C         STEGUN, AMS55, NATIONAL BUREAU OF STANDARDS, 1965, CHAPTER 9.
C
C         ASYMPTOTICS AND SPECIAL FUNCTIONS BY F.W.J. OLVER, ACADEMIC
C         PRESS, N.Y., 1974, PAGE 420
C
C     ABSTRACT
C         ZUNHJ COMPUTES PARAMETERS FOR BESSEL FUNCTIONS C(FNU,Z) =
C         J(FNU,Z), Y(FNU,Z) OR H(I,FNU,Z) I=1,2 FOR LARGE ORDERS FNU
C         BY MEANS OF THE UNIFORM ASYMPTOTIC EXPANSION
C
C         C(FNU,Z)=C1*PHI*( ASUM*AIRY(ARG) + C2*BSUM*DAIRY(ARG) )
C
C         FOR PROPER CHOICES OF C1, C2, AIRY AND DAIRY WHERE AIRY IS
C         AN AIRY FUNCTION AND DAIRY IS ITS DERIVATIVE.
C
C               (2/3)*FNU*ZETA**1.5 = ZETA1-ZETA2,
C
C         ZETA1=0.5*FNU*CLOG((1+W)/(1-W)), ZETA2=FNU*W FOR SCALING
C         PURPOSES IN AIRY FUNCTIONS FROM CAIRY OR CBIRY.
C
C         MCONJ=SIGN OF AIMAG(Z), BUT IS AMBIGUOUS WHEN Z IS REAL AND
C         MUST BE SPECIFIED. IPMTR=0 RETURNS ALL PARAMETERS. IPMTR=
C         1 COMPUTES ALL EXCEPT ASUM AND BSUM.
C
C***SEE ALSO  ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZDIV, ZLOG, ZSQRT
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   930122  Added ZLOG and ZSQRT to EXTERNAL statement.  (RWC)
C***END PROLOGUE  ZUNHJ
C     COMPLEX ARG,ASUM,BSUM,CFNU,CONE,CR,CZERO,DR,P,PHI,PRZTH,PTFN,
C    *RFN13,RTZTA,RZTH,SUMA,SUMB,TFN,T2,UP,W,W2,Z,ZA,ZB,ZC,ZETA,ZETA1,
C    *ZETA2,ZTH
      DOUBLE PRECISION ALFA, ANG, AP, AR, ARGI, ARGR, ASUMI, ASUMR,
     * ATOL, AW2, AZTH, BETA, BR, BSUMI, BSUMR, BTOL, C, CONEI, CONER,
     * CRI, CRR, DRI, DRR, EX1, EX2, FNU, FN13, FN23, GAMA, GPI, HPI,
     * PHII, PHIR, PI, PP, PR, PRZTHI, PRZTHR, PTFNI, PTFNR, RAW, RAW2,
     * RAZTH, RFNU, RFNU2, RFN13, RTZTI, RTZTR, RZTHI, RZTHR, STI, STR,
     * SUMAI, SUMAR, SUMBI, SUMBR, TEST, TFNI, TFNR, THPI, TOL, TZAI,
     * TZAR, T2I, T2R, UPI, UPR, WI, WR, W2I, W2R, ZAI, ZAR, ZBI, ZBR,
     * ZCI, ZCR, ZEROI, ZEROR, ZETAI, ZETAR, ZETA1I, ZETA1R, ZETA2I,
     * ZETA2R, ZI, ZR, ZTHI, ZTHR, ZABS, AC, D1MACH
      INTEGER IAS, IBS, IPMTR, IS, J, JR, JU, K, KMAX, KP1, KS, L, LR,
     * LRP1, L1, L2, M, IDUM
      DIMENSION AR(14), BR(14), C(105), ALFA(180), BETA(210), GAMA(30),
     * AP(30), PR(30), PI(30), UPR(14), UPI(14), CRR(14), CRI(14),
     * DRR(14), DRI(14)
      EXTERNAL ZABS, ZLOG, ZSQRT
      DATA AR(1), AR(2), AR(3), AR(4), AR(5), AR(6), AR(7), AR(8),
     1     AR(9), AR(10), AR(11), AR(12), AR(13), AR(14)/
     2     1.00000000000000000D+00,     1.04166666666666667D-01,
     3     8.35503472222222222D-02,     1.28226574556327160D-01,
     4     2.91849026464140464D-01,     8.81627267443757652D-01,
     5     3.32140828186276754D+00,     1.49957629868625547D+01,
     6     7.89230130115865181D+01,     4.74451538868264323D+02,
     7     3.20749009089066193D+03,     2.40865496408740049D+04,
     8     1.98923119169509794D+05,     1.79190200777534383D+06/
      DATA BR(1), BR(2), BR(3), BR(4), BR(5), BR(6), BR(7), BR(8),
     1     BR(9), BR(10), BR(11), BR(12), BR(13), BR(14)/
     2     1.00000000000000000D+00,    -1.45833333333333333D-01,
     3    -9.87413194444444444D-02,    -1.43312053915895062D-01,
     4    -3.17227202678413548D-01,    -9.42429147957120249D-01,
     5    -3.51120304082635426D+00,    -1.57272636203680451D+01,
     6    -8.22814390971859444D+01,    -4.92355370523670524D+02,
     7    -3.31621856854797251D+03,    -2.48276742452085896D+04,
     8    -2.04526587315129788D+05,    -1.83844491706820990D+06/
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2     C(19), C(20), C(21), C(22), C(23), C(24)/
     3     1.00000000000000000D+00,    -2.08333333333333333D-01,
     4     1.25000000000000000D-01,     3.34201388888888889D-01,
     5    -4.01041666666666667D-01,     7.03125000000000000D-02,
     6    -1.02581259645061728D+00,     1.84646267361111111D+00,
     7    -8.91210937500000000D-01,     7.32421875000000000D-02,
     8     4.66958442342624743D+00,    -1.12070026162229938D+01,
     9     8.78912353515625000D+00,    -2.36408691406250000D+00,
     A     1.12152099609375000D-01,    -2.82120725582002449D+01,
     B     8.46362176746007346D+01,    -9.18182415432400174D+01,
     C     4.25349987453884549D+01,    -7.36879435947963170D+00,
     D     2.27108001708984375D-01,     2.12570130039217123D+02,
     E    -7.65252468141181642D+02,     1.05999045252799988D+03/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3    -6.99579627376132541D+02,     2.18190511744211590D+02,
     4    -2.64914304869515555D+01,     5.72501420974731445D-01,
     5    -1.91945766231840700D+03,     8.06172218173730938D+03,
     6    -1.35865500064341374D+04,     1.16553933368645332D+04,
     7    -5.30564697861340311D+03,     1.20090291321635246D+03,
     8    -1.08090919788394656D+02,     1.72772750258445740D+00,
     9     2.02042913309661486D+04,    -9.69805983886375135D+04,
     A     1.92547001232531532D+05,    -2.03400177280415534D+05,
     B     1.22200464983017460D+05,    -4.11926549688975513D+04,
     C     7.10951430248936372D+03,    -4.93915304773088012D+02,
     D     6.07404200127348304D+00,    -2.42919187900551333D+05,
     E     1.31176361466297720D+06,    -2.99801591853810675D+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2     C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/
     3     3.76327129765640400D+06,    -2.81356322658653411D+06,
     4     1.26836527332162478D+06,    -3.31645172484563578D+05,
     5     4.52187689813627263D+04,    -2.49983048181120962D+03,
     6     2.43805296995560639D+01,     3.28446985307203782D+06,
     7    -1.97068191184322269D+07,     5.09526024926646422D+07,
     8    -7.41051482115326577D+07,     6.63445122747290267D+07,
     9    -3.75671766607633513D+07,     1.32887671664218183D+07,
     A    -2.78561812808645469D+06,     3.08186404612662398D+05,
     B    -1.38860897537170405D+04,     1.10017140269246738D+02,
     C    -4.93292536645099620D+07,     3.25573074185765749D+08,
     D    -9.39462359681578403D+08,     1.55359689957058006D+09,
     E    -1.62108055210833708D+09,     1.10684281682301447D+09/
      DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),
     1     C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),
     2     C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/
     3    -4.95889784275030309D+08,     1.42062907797533095D+08,
     4    -2.44740627257387285D+07,     2.24376817792244943D+06,
     5    -8.40054336030240853D+04,     5.51335896122020586D+02,
     6     8.14789096118312115D+08,    -5.86648149205184723D+09,
     7     1.86882075092958249D+10,    -3.46320433881587779D+10,
     8     4.12801855797539740D+10,    -3.30265997498007231D+10,
     9     1.79542137311556001D+10,    -6.56329379261928433D+09,
     A     1.55927986487925751D+09,    -2.25105661889415278D+08,
     B     1.73951075539781645D+07,    -5.49842327572288687D+05,
     C     3.03809051092238427D+03,    -1.46792612476956167D+10,
     D     1.14498237732025810D+11,    -3.99096175224466498D+11,
     E     8.19218669548577329D+11,    -1.09837515608122331D+12/
      DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104),
     1     C(105)/
     2     1.00815810686538209D+12,    -6.45364869245376503D+11,
     3     2.87900649906150589D+11,    -8.78670721780232657D+10,
     4     1.76347306068349694D+10,    -2.16716498322379509D+09,
     5     1.43157876718888981D+08,    -3.87183344257261262D+06,
     6     1.82577554742931747D+04/
      DATA ALFA(1), ALFA(2), ALFA(3), ALFA(4), ALFA(5), ALFA(6),
     1     ALFA(7), ALFA(8), ALFA(9), ALFA(10), ALFA(11), ALFA(12),
     2     ALFA(13), ALFA(14), ALFA(15), ALFA(16), ALFA(17), ALFA(18),
     3     ALFA(19), ALFA(20), ALFA(21), ALFA(22)/
     4    -4.44444444444444444D-03,    -9.22077922077922078D-04,
     5    -8.84892884892884893D-05,     1.65927687832449737D-04,
     6     2.46691372741792910D-04,     2.65995589346254780D-04,
     7     2.61824297061500945D-04,     2.48730437344655609D-04,
     8     2.32721040083232098D-04,     2.16362485712365082D-04,
     9     2.00738858762752355D-04,     1.86267636637545172D-04,
     A     1.73060775917876493D-04,     1.61091705929015752D-04,
     B     1.50274774160908134D-04,     1.40503497391269794D-04,
     C     1.31668816545922806D-04,     1.23667445598253261D-04,
     D     1.16405271474737902D-04,     1.09798298372713369D-04,
     E     1.03772410422992823D-04,     9.82626078369363448D-05/
      DATA ALFA(23), ALFA(24), ALFA(25), ALFA(26), ALFA(27), ALFA(28),
     1     ALFA(29), ALFA(30), ALFA(31), ALFA(32), ALFA(33), ALFA(34),
     2     ALFA(35), ALFA(36), ALFA(37), ALFA(38), ALFA(39), ALFA(40),
     3     ALFA(41), ALFA(42), ALFA(43), ALFA(44)/
     4     9.32120517249503256D-05,     8.85710852478711718D-05,
     5     8.42963105715700223D-05,     8.03497548407791151D-05,
     6     7.66981345359207388D-05,     7.33122157481777809D-05,
     7     7.01662625163141333D-05,     6.72375633790160292D-05,
     8     6.93735541354588974D-04,     2.32241745182921654D-04,
     9    -1.41986273556691197D-05,    -1.16444931672048640D-04,
     A    -1.50803558053048762D-04,    -1.55121924918096223D-04,
     B    -1.46809756646465549D-04,    -1.33815503867491367D-04,
     C    -1.19744975684254051D-04,    -1.06184319207974020D-04,
     D    -9.37699549891194492D-05,    -8.26923045588193274D-05,
     E    -7.29374348155221211D-05,    -6.44042357721016283D-05/
      DATA ALFA(45), ALFA(46), ALFA(47), ALFA(48), ALFA(49), ALFA(50),
     1     ALFA(51), ALFA(52), ALFA(53), ALFA(54), ALFA(55), ALFA(56),
     2     ALFA(57), ALFA(58), ALFA(59), ALFA(60), ALFA(61), ALFA(62),
     3     ALFA(63), ALFA(64), ALFA(65), ALFA(66)/
     4    -5.69611566009369048D-05,    -5.04731044303561628D-05,
     5    -4.48134868008882786D-05,    -3.98688727717598864D-05,
     6    -3.55400532972042498D-05,    -3.17414256609022480D-05,
     7    -2.83996793904174811D-05,    -2.54522720634870566D-05,
     8    -2.28459297164724555D-05,    -2.05352753106480604D-05,
     9    -1.84816217627666085D-05,    -1.66519330021393806D-05,
     A    -1.50179412980119482D-05,    -1.35554031379040526D-05,
     B    -1.22434746473858131D-05,    -1.10641884811308169D-05,
     C    -3.54211971457743841D-04,    -1.56161263945159416D-04,
     D     3.04465503594936410D-05,     1.30198655773242693D-04,
     E     1.67471106699712269D-04,     1.70222587683592569D-04/
      DATA ALFA(67), ALFA(68), ALFA(69), ALFA(70), ALFA(71), ALFA(72),
     1     ALFA(73), ALFA(74), ALFA(75), ALFA(76), ALFA(77), ALFA(78),
     2     ALFA(79), ALFA(80), ALFA(81), ALFA(82), ALFA(83), ALFA(84),
     3     ALFA(85), ALFA(86), ALFA(87), ALFA(88)/
     4     1.56501427608594704D-04,     1.36339170977445120D-04,
     5     1.14886692029825128D-04,     9.45869093034688111D-05,
     6     7.64498419250898258D-05,     6.07570334965197354D-05,
     7     4.74394299290508799D-05,     3.62757512005344297D-05,
     8     2.69939714979224901D-05,     1.93210938247939253D-05,
     9     1.30056674793963203D-05,     7.82620866744496661D-06,
     A     3.59257485819351583D-06,     1.44040049814251817D-07,
     B    -2.65396769697939116D-06,    -4.91346867098485910D-06,
     C    -6.72739296091248287D-06,    -8.17269379678657923D-06,
     D    -9.31304715093561232D-06,    -1.02011418798016441D-05,
     E    -1.08805962510592880D-05,    -1.13875481509603555D-05/
      DATA ALFA(89), ALFA(90), ALFA(91), ALFA(92), ALFA(93), ALFA(94),
     1     ALFA(95), ALFA(96), ALFA(97), ALFA(98), ALFA(99), ALFA(100),
     2     ALFA(101), ALFA(102), ALFA(103), ALFA(104), ALFA(105),
     3     ALFA(106), ALFA(107), ALFA(108), ALFA(109), ALFA(110)/
     4    -1.17519675674556414D-05,    -1.19987364870944141D-05,
     5     3.78194199201772914D-04,     2.02471952761816167D-04,
     6    -6.37938506318862408D-05,    -2.38598230603005903D-04,
     7    -3.10916256027361568D-04,    -3.13680115247576316D-04,
     8    -2.78950273791323387D-04,    -2.28564082619141374D-04,
     9    -1.75245280340846749D-04,    -1.25544063060690348D-04,
     A    -8.22982872820208365D-05,    -4.62860730588116458D-05,
     B    -1.72334302366962267D-05,     5.60690482304602267D-06,
     C     2.31395443148286800D-05,     3.62642745856793957D-05,
     D     4.58006124490188752D-05,     5.24595294959114050D-05,
     E     5.68396208545815266D-05,     5.94349820393104052D-05/
      DATA ALFA(111), ALFA(112), ALFA(113), ALFA(114), ALFA(115),
     1     ALFA(116), ALFA(117), ALFA(118), ALFA(119), ALFA(120),
     2     ALFA(121), ALFA(122), ALFA(123), ALFA(124), ALFA(125),
     3     ALFA(126), ALFA(127), ALFA(128), ALFA(129), ALFA(130)/
     4     6.06478527578421742D-05,     6.08023907788436497D-05,
     5     6.01577894539460388D-05,     5.89199657344698500D-05,
     6     5.72515823777593053D-05,     5.52804375585852577D-05,
     7     5.31063773802880170D-05,     5.08069302012325706D-05,
     8     4.84418647620094842D-05,     4.60568581607475370D-05,
     9    -6.91141397288294174D-04,    -4.29976633058871912D-04,
     A     1.83067735980039018D-04,     6.60088147542014144D-04,
     B     8.75964969951185931D-04,     8.77335235958235514D-04,
     C     7.49369585378990637D-04,     5.63832329756980918D-04,
     D     3.68059319971443156D-04,     1.88464535514455599D-04/
      DATA ALFA(131), ALFA(132), ALFA(133), ALFA(134), ALFA(135),
     1     ALFA(136), ALFA(137), ALFA(138), ALFA(139), ALFA(140),
     2     ALFA(141), ALFA(142), ALFA(143), ALFA(144), ALFA(145),
     3     ALFA(146), ALFA(147), ALFA(148), ALFA(149), ALFA(150)/
     4     3.70663057664904149D-05,    -8.28520220232137023D-05,
     5    -1.72751952869172998D-04,    -2.36314873605872983D-04,
     6    -2.77966150694906658D-04,    -3.02079514155456919D-04,
     7    -3.12594712643820127D-04,    -3.12872558758067163D-04,
     8    -3.05678038466324377D-04,    -2.93226470614557331D-04,
     9    -2.77255655582934777D-04,    -2.59103928467031709D-04,
     A    -2.39784014396480342D-04,    -2.20048260045422848D-04,
     B    -2.00443911094971498D-04,    -1.81358692210970687D-04,
     C    -1.63057674478657464D-04,    -1.45712672175205844D-04,
     D    -1.29425421983924587D-04,    -1.14245691942445952D-04/
      DATA ALFA(151), ALFA(152), ALFA(153), ALFA(154), ALFA(155),
     1     ALFA(156), ALFA(157), ALFA(158), ALFA(159), ALFA(160),
     2     ALFA(161), ALFA(162), ALFA(163), ALFA(164), ALFA(165),
     3     ALFA(166), ALFA(167), ALFA(168), ALFA(169), ALFA(170)/
     4     1.92821964248775885D-03,     1.35592576302022234D-03,
     5    -7.17858090421302995D-04,    -2.58084802575270346D-03,
     6    -3.49271130826168475D-03,    -3.46986299340960628D-03,
     7    -2.82285233351310182D-03,    -1.88103076404891354D-03,
     8    -8.89531718383947600D-04,     3.87912102631035228D-06,
     9     7.28688540119691412D-04,     1.26566373053457758D-03,
     A     1.62518158372674427D-03,     1.83203153216373172D-03,
     B     1.91588388990527909D-03,     1.90588846755546138D-03,
     C     1.82798982421825727D-03,     1.70389506421121530D-03,
     D     1.55097127171097686D-03,     1.38261421852276159D-03/
      DATA ALFA(171), ALFA(172), ALFA(173), ALFA(174), ALFA(175),
     1     ALFA(176), ALFA(177), ALFA(178), ALFA(179), ALFA(180)/
     2     1.20881424230064774D-03,     1.03676532638344962D-03,
     3     8.71437918068619115D-04,     7.16080155297701002D-04,
     4     5.72637002558129372D-04,     4.42089819465802277D-04,
     5     3.24724948503090564D-04,     2.20342042730246599D-04,
     6     1.28412898401353882D-04,     4.82005924552095464D-05/
      DATA BETA(1), BETA(2), BETA(3), BETA(4), BETA(5), BETA(6),
     1     BETA(7), BETA(8), BETA(9), BETA(10), BETA(11), BETA(12),
     2     BETA(13), BETA(14), BETA(15), BETA(16), BETA(17), BETA(18),
     3     BETA(19), BETA(20), BETA(21), BETA(22)/
     4     1.79988721413553309D-02,     5.59964911064388073D-03,
     5     2.88501402231132779D-03,     1.80096606761053941D-03,
     6     1.24753110589199202D-03,     9.22878876572938311D-04,
     7     7.14430421727287357D-04,     5.71787281789704872D-04,
     8     4.69431007606481533D-04,     3.93232835462916638D-04,
     9     3.34818889318297664D-04,     2.88952148495751517D-04,
     A     2.52211615549573284D-04,     2.22280580798883327D-04,
     B     1.97541838033062524D-04,     1.76836855019718004D-04,
     C     1.59316899661821081D-04,     1.44347930197333986D-04,
     D     1.31448068119965379D-04,     1.20245444949302884D-04,
     E     1.10449144504599392D-04,     1.01828770740567258D-04/
      DATA BETA(23), BETA(24), BETA(25), BETA(26), BETA(27), BETA(28),
     1     BETA(29), BETA(30), BETA(31), BETA(32), BETA(33), BETA(34),
     2     BETA(35), BETA(36), BETA(37), BETA(38), BETA(39), BETA(40),
     3     BETA(41), BETA(42), BETA(43), BETA(44)/
     4     9.41998224204237509D-05,     8.74130545753834437D-05,
     5     8.13466262162801467D-05,     7.59002269646219339D-05,
     6     7.09906300634153481D-05,     6.65482874842468183D-05,
     7     6.25146958969275078D-05,     5.88403394426251749D-05,
     8    -1.49282953213429172D-03,    -8.78204709546389328D-04,
     9    -5.02916549572034614D-04,    -2.94822138512746025D-04,
     A    -1.75463996970782828D-04,    -1.04008550460816434D-04,
     B    -5.96141953046457895D-05,    -3.12038929076098340D-05,
     C    -1.26089735980230047D-05,    -2.42892608575730389D-07,
     D     8.05996165414273571D-06,     1.36507009262147391D-05,
     E     1.73964125472926261D-05,     1.98672978842133780D-05/
      DATA BETA(45), BETA(46), BETA(47), BETA(48), BETA(49), BETA(50),
     1     BETA(51), BETA(52), BETA(53), BETA(54), BETA(55), BETA(56),
     2     BETA(57), BETA(58), BETA(59), BETA(60), BETA(61), BETA(62),
     3     BETA(63), BETA(64), BETA(65), BETA(66)/
     4     2.14463263790822639D-05,     2.23954659232456514D-05,
     5     2.28967783814712629D-05,     2.30785389811177817D-05,
     6     2.30321976080909144D-05,     2.28236073720348722D-05,
     7     2.25005881105292418D-05,     2.20981015361991429D-05,
     8     2.16418427448103905D-05,     2.11507649256220843D-05,
     9     2.06388749782170737D-05,     2.01165241997081666D-05,
     A     1.95913450141179244D-05,     1.90689367910436740D-05,
     B     1.85533719641636667D-05,     1.80475722259674218D-05,
     C     5.52213076721292790D-04,     4.47932581552384646D-04,
     D     2.79520653992020589D-04,     1.52468156198446602D-04,
     E     6.93271105657043598D-05,     1.76258683069991397D-05/
      DATA BETA(67), BETA(68), BETA(69), BETA(70), BETA(71), BETA(72),
     1     BETA(73), BETA(74), BETA(75), BETA(76), BETA(77), BETA(78),
     2     BETA(79), BETA(80), BETA(81), BETA(82), BETA(83), BETA(84),
     3     BETA(85), BETA(86), BETA(87), BETA(88)/
     4    -1.35744996343269136D-05,    -3.17972413350427135D-05,
     5    -4.18861861696693365D-05,    -4.69004889379141029D-05,
     6    -4.87665447413787352D-05,    -4.87010031186735069D-05,
     7    -4.74755620890086638D-05,    -4.55813058138628452D-05,
     8    -4.33309644511266036D-05,    -4.09230193157750364D-05,
     9    -3.84822638603221274D-05,    -3.60857167535410501D-05,
     A    -3.37793306123367417D-05,    -3.15888560772109621D-05,
     B    -2.95269561750807315D-05,    -2.75978914828335759D-05,
     C    -2.58006174666883713D-05,    -2.41308356761280200D-05,
     D    -2.25823509518346033D-05,    -2.11479656768912971D-05,
     E    -1.98200638885294927D-05,    -1.85909870801065077D-05/
      DATA BETA(89), BETA(90), BETA(91), BETA(92), BETA(93), BETA(94),
     1     BETA(95), BETA(96), BETA(97), BETA(98), BETA(99), BETA(100),
     2     BETA(101), BETA(102), BETA(103), BETA(104), BETA(105),
     3     BETA(106), BETA(107), BETA(108), BETA(109), BETA(110)/
     4    -1.74532699844210224D-05,    -1.63997823854497997D-05,
     5    -4.74617796559959808D-04,    -4.77864567147321487D-04,
     6    -3.20390228067037603D-04,    -1.61105016119962282D-04,
     7    -4.25778101285435204D-05,     3.44571294294967503D-05,
     8     7.97092684075674924D-05,     1.03138236708272200D-04,
     9     1.12466775262204158D-04,     1.13103642108481389D-04,
     A     1.08651634848774268D-04,     1.01437951597661973D-04,
     B     9.29298396593363896D-05,     8.40293133016089978D-05,
     C     7.52727991349134062D-05,     6.69632521975730872D-05,
     D     5.92564547323194704D-05,     5.22169308826975567D-05,
     E     4.58539485165360646D-05,     4.01445513891486808D-05/
      DATA BETA(111), BETA(112), BETA(113), BETA(114), BETA(115),
     1     BETA(116), BETA(117), BETA(118), BETA(119), BETA(120),
     2     BETA(121), BETA(122), BETA(123), BETA(124), BETA(125),
     3     BETA(126), BETA(127), BETA(128), BETA(129), BETA(130)/
     4     3.50481730031328081D-05,     3.05157995034346659D-05,
     5     2.64956119950516039D-05,     2.29363633690998152D-05,
     6     1.97893056664021636D-05,     1.70091984636412623D-05,
     7     1.45547428261524004D-05,     1.23886640995878413D-05,
     8     1.04775876076583236D-05,     8.79179954978479373D-06,
     9     7.36465810572578444D-04,     8.72790805146193976D-04,
     A     6.22614862573135066D-04,     2.85998154194304147D-04,
     B     3.84737672879366102D-06,    -1.87906003636971558D-04,
     C    -2.97603646594554535D-04,    -3.45998126832656348D-04,
     D    -3.53382470916037712D-04,    -3.35715635775048757D-04/
      DATA BETA(131), BETA(132), BETA(133), BETA(134), BETA(135),
     1     BETA(136), BETA(137), BETA(138), BETA(139), BETA(140),
     2     BETA(141), BETA(142), BETA(143), BETA(144), BETA(145),
     3     BETA(146), BETA(147), BETA(148), BETA(149), BETA(150)/
     4    -3.04321124789039809D-04,    -2.66722723047612821D-04,
     5    -2.27654214122819527D-04,    -1.89922611854562356D-04,
     6    -1.55058918599093870D-04,    -1.23778240761873630D-04,
     7    -9.62926147717644187D-05,    -7.25178327714425337D-05,
     8    -5.22070028895633801D-05,    -3.50347750511900522D-05,
     9    -2.06489761035551757D-05,    -8.70106096849767054D-06,
     A     1.13698686675100290D-06,     9.16426474122778849D-06,
     B     1.56477785428872620D-05,     2.08223629482466847D-05,
     C     2.48923381004595156D-05,     2.80340509574146325D-05,
     D     3.03987774629861915D-05,     3.21156731406700616D-05/
      DATA BETA(151), BETA(152), BETA(153), BETA(154), BETA(155),
     1     BETA(156), BETA(157), BETA(158), BETA(159), BETA(160),
     2     BETA(161), BETA(162), BETA(163), BETA(164), BETA(165),
     3     BETA(166), BETA(167), BETA(168), BETA(169), BETA(170)/
     4    -1.80182191963885708D-03,    -2.43402962938042533D-03,
     5    -1.83422663549856802D-03,    -7.62204596354009765D-04,
     6     2.39079475256927218D-04,     9.49266117176881141D-04,
     7     1.34467449701540359D-03,     1.48457495259449178D-03,
     8     1.44732339830617591D-03,     1.30268261285657186D-03,
     9     1.10351597375642682D-03,     8.86047440419791759D-04,
     A     6.73073208165665473D-04,     4.77603872856582378D-04,
     B     3.05991926358789362D-04,     1.60315694594721630D-04,
     C     4.00749555270613286D-05,    -5.66607461635251611D-05,
     D    -1.32506186772982638D-04,    -1.90296187989614057D-04/
      DATA BETA(171), BETA(172), BETA(173), BETA(174), BETA(175),
     1     BETA(176), BETA(177), BETA(178), BETA(179), BETA(180),
     2     BETA(181), BETA(182), BETA(183), BETA(184), BETA(185),
     3     BETA(186), BETA(187), BETA(188), BETA(189), BETA(190)/
     4    -2.32811450376937408D-04,    -2.62628811464668841D-04,
     5    -2.82050469867598672D-04,    -2.93081563192861167D-04,
     6    -2.97435962176316616D-04,    -2.96557334239348078D-04,
     7    -2.91647363312090861D-04,    -2.83696203837734166D-04,
     8    -2.73512317095673346D-04,    -2.61750155806768580D-04,
     9     6.38585891212050914D-03,     9.62374215806377941D-03,
     A     7.61878061207001043D-03,     2.83219055545628054D-03,
     B    -2.09841352012720090D-03,    -5.73826764216626498D-03,
     C    -7.70804244495414620D-03,    -8.21011692264844401D-03,
     D    -7.65824520346905413D-03,    -6.47209729391045177D-03/
      DATA BETA(191), BETA(192), BETA(193), BETA(194), BETA(195),
     1     BETA(196), BETA(197), BETA(198), BETA(199), BETA(200),
     2     BETA(201), BETA(202), BETA(203), BETA(204), BETA(205),
     3     BETA(206), BETA(207), BETA(208), BETA(209), BETA(210)/
     4    -4.99132412004966473D-03,    -3.45612289713133280D-03,
     5    -2.01785580014170775D-03,    -7.59430686781961401D-04,
     6     2.84173631523859138D-04,     1.10891667586337403D-03,
     7     1.72901493872728771D-03,     2.16812590802684701D-03,
     8     2.45357710494539735D-03,     2.61281821058334862D-03,
     9     2.67141039656276912D-03,     2.65203073395980430D-03,
     A     2.57411652877287315D-03,     2.45389126236094427D-03,
     B     2.30460058071795494D-03,     2.13684837686712662D-03,
     C     1.95896528478870911D-03,     1.77737008679454412D-03,
     D     1.59690280765839059D-03,     1.42111975664438546D-03/
      DATA GAMA(1), GAMA(2), GAMA(3), GAMA(4), GAMA(5), GAMA(6),
     1     GAMA(7), GAMA(8), GAMA(9), GAMA(10), GAMA(11), GAMA(12),
     2     GAMA(13), GAMA(14), GAMA(15), GAMA(16), GAMA(17), GAMA(18),
     3     GAMA(19), GAMA(20), GAMA(21), GAMA(22)/
     4     6.29960524947436582D-01,     2.51984209978974633D-01,
     5     1.54790300415655846D-01,     1.10713062416159013D-01,
     6     8.57309395527394825D-02,     6.97161316958684292D-02,
     7     5.86085671893713576D-02,     5.04698873536310685D-02,
     8     4.42600580689154809D-02,     3.93720661543509966D-02,
     9     3.54283195924455368D-02,     3.21818857502098231D-02,
     A     2.94646240791157679D-02,     2.71581677112934479D-02,
     B     2.51768272973861779D-02,     2.34570755306078891D-02,
     C     2.19508390134907203D-02,     2.06210828235646240D-02,
     D     1.94388240897880846D-02,     1.83810633800683158D-02,
     E     1.74293213231963172D-02,     1.65685837786612353D-02/
      DATA GAMA(23), GAMA(24), GAMA(25), GAMA(26), GAMA(27), GAMA(28),
     1     GAMA(29), GAMA(30)/
     2     1.57865285987918445D-02,     1.50729501494095594D-02,
     3     1.44193250839954639D-02,     1.38184805735341786D-02,
     4     1.32643378994276568D-02,     1.27517121970498651D-02,
     5     1.22761545318762767D-02,     1.18338262398482403D-02/
      DATA EX1, EX2, HPI, GPI, THPI /
     1     3.33333333333333333D-01,     6.66666666666666667D-01,
     2     1.57079632679489662D+00,     3.14159265358979324D+00,
     3     4.71238898038468986D+00/
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZUNHJ
      RFNU = 1.0D0/FNU
C-----------------------------------------------------------------------
C     OVERFLOW TEST (Z/FNU TOO SMALL)
C-----------------------------------------------------------------------
      TEST = D1MACH(1)*1.0D+3
      AC = FNU*TEST
      IF (ABS(ZR).GT.AC .OR. ABS(ZI).GT.AC) GO TO 15
      ZETA1R = 2.0D0*ABS(LOG(TEST))+FNU
      ZETA1I = 0.0D0
      ZETA2R = FNU
      ZETA2I = 0.0D0
      PHIR = 1.0D0
      PHII = 0.0D0
      ARGR = 1.0D0
      ARGI = 0.0D0
      RETURN
   15 CONTINUE
      ZBR = ZR*RFNU
      ZBI = ZI*RFNU
      RFNU2 = RFNU*RFNU
C-----------------------------------------------------------------------
C     COMPUTE IN THE FOURTH QUADRANT
C-----------------------------------------------------------------------
      FN13 = FNU**EX1
      FN23 = FN13*FN13
      RFN13 = 1.0D0/FN13
      W2R = CONER - ZBR*ZBR + ZBI*ZBI
      W2I = CONEI - ZBR*ZBI - ZBR*ZBI
      AW2 = ZABS(W2R,W2I)
      IF (AW2.GT.0.25D0) GO TO 130
C-----------------------------------------------------------------------
C     POWER SERIES FOR ABS(W2).LE.0.25D0
C-----------------------------------------------------------------------
      K = 1
      PR(1) = CONER
      PI(1) = CONEI
      SUMAR = GAMA(1)
      SUMAI = ZEROI
      AP(1) = 1.0D0
      IF (AW2.LT.TOL) GO TO 20
      DO 10 K=2,30
        PR(K) = PR(K-1)*W2R - PI(K-1)*W2I
        PI(K) = PR(K-1)*W2I + PI(K-1)*W2R
        SUMAR = SUMAR + PR(K)*GAMA(K)
        SUMAI = SUMAI + PI(K)*GAMA(K)
        AP(K) = AP(K-1)*AW2
        IF (AP(K).LT.TOL) GO TO 20
   10 CONTINUE
      K = 30
   20 CONTINUE
      KMAX = K
      ZETAR = W2R*SUMAR - W2I*SUMAI
      ZETAI = W2R*SUMAI + W2I*SUMAR
      ARGR = ZETAR*FN23
      ARGI = ZETAI*FN23
      CALL ZSQRT(SUMAR, SUMAI, ZAR, ZAI)
      CALL ZSQRT(W2R, W2I, STR, STI)
      ZETA2R = STR*FNU
      ZETA2I = STI*FNU
      STR = CONER + EX2*(ZETAR*ZAR-ZETAI*ZAI)
      STI = CONEI + EX2*(ZETAR*ZAI+ZETAI*ZAR)
      ZETA1R = STR*ZETA2R - STI*ZETA2I
      ZETA1I = STR*ZETA2I + STI*ZETA2R
      ZAR = ZAR + ZAR
      ZAI = ZAI + ZAI
      CALL ZSQRT(ZAR, ZAI, STR, STI)
      PHIR = STR*RFN13
      PHII = STI*RFN13
      IF (IPMTR.EQ.1) GO TO 120
C-----------------------------------------------------------------------
C     SUM SERIES FOR ASUM AND BSUM
C-----------------------------------------------------------------------
      SUMBR = ZEROR
      SUMBI = ZEROI
      DO 30 K=1,KMAX
        SUMBR = SUMBR + PR(K)*BETA(K)
        SUMBI = SUMBI + PI(K)*BETA(K)
   30 CONTINUE
      ASUMR = ZEROR
      ASUMI = ZEROI
      BSUMR = SUMBR
      BSUMI = SUMBI
      L1 = 0
      L2 = 30
      BTOL = TOL*(ABS(BSUMR)+ABS(BSUMI))
      ATOL = TOL
      PP = 1.0D0
      IAS = 0
      IBS = 0
      IF (RFNU2.LT.TOL) GO TO 110
      DO 100 IS=2,7
        ATOL = ATOL/RFNU2
        PP = PP*RFNU2
        IF (IAS.EQ.1) GO TO 60
        SUMAR = ZEROR
        SUMAI = ZEROI
        DO 40 K=1,KMAX
          M = L1 + K
          SUMAR = SUMAR + PR(K)*ALFA(M)
          SUMAI = SUMAI + PI(K)*ALFA(M)
          IF (AP(K).LT.ATOL) GO TO 50
   40   CONTINUE
   50   CONTINUE
        ASUMR = ASUMR + SUMAR*PP
        ASUMI = ASUMI + SUMAI*PP
        IF (PP.LT.TOL) IAS = 1
   60   CONTINUE
        IF (IBS.EQ.1) GO TO 90
        SUMBR = ZEROR
        SUMBI = ZEROI
        DO 70 K=1,KMAX
          M = L2 + K
          SUMBR = SUMBR + PR(K)*BETA(M)
          SUMBI = SUMBI + PI(K)*BETA(M)
          IF (AP(K).LT.ATOL) GO TO 80
   70   CONTINUE
   80   CONTINUE
        BSUMR = BSUMR + SUMBR*PP
        BSUMI = BSUMI + SUMBI*PP
        IF (PP.LT.BTOL) IBS = 1
   90   CONTINUE
        IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 110
        L1 = L1 + 30
        L2 = L2 + 30
  100 CONTINUE
  110 CONTINUE
      ASUMR = ASUMR + CONER
      PP = RFNU*RFN13
      BSUMR = BSUMR*PP
      BSUMI = BSUMI*PP
  120 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     ABS(W2).GT.0.25D0
C-----------------------------------------------------------------------
  130 CONTINUE
      CALL ZSQRT(W2R, W2I, WR, WI)
      IF (WR.LT.0.0D0) WR = 0.0D0
      IF (WI.LT.0.0D0) WI = 0.0D0
      STR = CONER + WR
      STI = WI
      CALL ZDIV(STR, STI, ZBR, ZBI, ZAR, ZAI)
      CALL ZLOG(ZAR, ZAI, ZCR, ZCI, IDUM)
      IF (ZCI.LT.0.0D0) ZCI = 0.0D0
      IF (ZCI.GT.HPI) ZCI = HPI
      IF (ZCR.LT.0.0D0) ZCR = 0.0D0
      ZTHR = (ZCR-WR)*1.5D0
      ZTHI = (ZCI-WI)*1.5D0
      ZETA1R = ZCR*FNU
      ZETA1I = ZCI*FNU
      ZETA2R = WR*FNU
      ZETA2I = WI*FNU
      AZTH = ZABS(ZTHR,ZTHI)
      ANG = THPI
      IF (ZTHR.GE.0.0D0 .AND. ZTHI.LT.0.0D0) GO TO 140
      ANG = HPI
      IF (ZTHR.EQ.0.0D0) GO TO 140
      ANG = DATAN(ZTHI/ZTHR)
      IF (ZTHR.LT.0.0D0) ANG = ANG + GPI
  140 CONTINUE
      PP = AZTH**EX2
      ANG = ANG*EX2
      ZETAR = PP*COS(ANG)
      ZETAI = PP*SIN(ANG)
      IF (ZETAI.LT.0.0D0) ZETAI = 0.0D0
      ARGR = ZETAR*FN23
      ARGI = ZETAI*FN23
      CALL ZDIV(ZTHR, ZTHI, ZETAR, ZETAI, RTZTR, RTZTI)
      CALL ZDIV(RTZTR, RTZTI, WR, WI, ZAR, ZAI)
      TZAR = ZAR + ZAR
      TZAI = ZAI + ZAI
      CALL ZSQRT(TZAR, TZAI, STR, STI)
      PHIR = STR*RFN13
      PHII = STI*RFN13
      IF (IPMTR.EQ.1) GO TO 120
      RAW = 1.0D0/SQRT(AW2)
      STR = WR*RAW
      STI = -WI*RAW
      TFNR = STR*RFNU*RAW
      TFNI = STI*RFNU*RAW
      RAZTH = 1.0D0/AZTH
      STR = ZTHR*RAZTH
      STI = -ZTHI*RAZTH
      RZTHR = STR*RAZTH*RFNU
      RZTHI = STI*RAZTH*RFNU
      ZCR = RZTHR*AR(2)
      ZCI = RZTHI*AR(2)
      RAW2 = 1.0D0/AW2
      STR = W2R*RAW2
      STI = -W2I*RAW2
      T2R = STR*RAW2
      T2I = STI*RAW2
      STR = T2R*C(2) + C(3)
      STI = T2I*C(2)
      UPR(2) = STR*TFNR - STI*TFNI
      UPI(2) = STR*TFNI + STI*TFNR
      BSUMR = UPR(2) + ZCR
      BSUMI = UPI(2) + ZCI
      ASUMR = ZEROR
      ASUMI = ZEROI
      IF (RFNU.LT.TOL) GO TO 220
      PRZTHR = RZTHR
      PRZTHI = RZTHI
      PTFNR = TFNR
      PTFNI = TFNI
      UPR(1) = CONER
      UPI(1) = CONEI
      PP = 1.0D0
      BTOL = TOL*(ABS(BSUMR)+ABS(BSUMI))
      KS = 0
      KP1 = 2
      L = 3
      IAS = 0
      IBS = 0
      DO 210 LR=2,12,2
        LRP1 = LR + 1
C-----------------------------------------------------------------------
C     COMPUTE TWO ADDITIONAL CR, DR, AND UP FOR TWO MORE TERMS IN
C     NEXT SUMA AND SUMB
C-----------------------------------------------------------------------
        DO 160 K=LR,LRP1
          KS = KS + 1
          KP1 = KP1 + 1
          L = L + 1
          ZAR = C(L)
          ZAI = ZEROI
          DO 150 J=2,KP1
            L = L + 1
            STR = ZAR*T2R - T2I*ZAI + C(L)
            ZAI = ZAR*T2I + ZAI*T2R
            ZAR = STR
  150     CONTINUE
          STR = PTFNR*TFNR - PTFNI*TFNI
          PTFNI = PTFNR*TFNI + PTFNI*TFNR
          PTFNR = STR
          UPR(KP1) = PTFNR*ZAR - PTFNI*ZAI
          UPI(KP1) = PTFNI*ZAR + PTFNR*ZAI
          CRR(KS) = PRZTHR*BR(KS+1)
          CRI(KS) = PRZTHI*BR(KS+1)
          STR = PRZTHR*RZTHR - PRZTHI*RZTHI
          PRZTHI = PRZTHR*RZTHI + PRZTHI*RZTHR
          PRZTHR = STR
          DRR(KS) = PRZTHR*AR(KS+2)
          DRI(KS) = PRZTHI*AR(KS+2)
  160   CONTINUE
        PP = PP*RFNU2
        IF (IAS.EQ.1) GO TO 180
        SUMAR = UPR(LRP1)
        SUMAI = UPI(LRP1)
        JU = LRP1
        DO 170 JR=1,LR
          JU = JU - 1
          SUMAR = SUMAR + CRR(JR)*UPR(JU) - CRI(JR)*UPI(JU)
          SUMAI = SUMAI + CRR(JR)*UPI(JU) + CRI(JR)*UPR(JU)
  170   CONTINUE
        ASUMR = ASUMR + SUMAR
        ASUMI = ASUMI + SUMAI
        TEST = ABS(SUMAR) + ABS(SUMAI)
        IF (PP.LT.TOL .AND. TEST.LT.TOL) IAS = 1
  180   CONTINUE
        IF (IBS.EQ.1) GO TO 200
        SUMBR = UPR(LR+2) + UPR(LRP1)*ZCR - UPI(LRP1)*ZCI
        SUMBI = UPI(LR+2) + UPR(LRP1)*ZCI + UPI(LRP1)*ZCR
        JU = LRP1
        DO 190 JR=1,LR
          JU = JU - 1
          SUMBR = SUMBR + DRR(JR)*UPR(JU) - DRI(JR)*UPI(JU)
          SUMBI = SUMBI + DRR(JR)*UPI(JU) + DRI(JR)*UPR(JU)
  190   CONTINUE
        BSUMR = BSUMR + SUMBR
        BSUMI = BSUMI + SUMBI
        TEST = ABS(SUMBR) + ABS(SUMBI)
        IF (PP.LT.BTOL .AND. TEST.LT.BTOL) IBS = 1
  200   CONTINUE
        IF (IAS.EQ.1 .AND. IBS.EQ.1) GO TO 220
  210 CONTINUE
  220 CONTINUE
      ASUMR = ASUMR + CONER
      STR = -BSUMR*RFN13
      STI = -BSUMI*RFN13
      CALL ZDIV(STR, STI, RTZTR, RTZTI, BSUMR, BSUMI)
      GO TO 120
      END
*DECK ZUNI1
      SUBROUTINE ZUNI1 (ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL,
     +   TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  ZUNI1
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CUNI1-A, ZUNI1-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
C     EXPANSION FOR I(FNU,Z) IN -PI/3.LE.ARG Z.LE.PI/3.
C
C     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
C     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
C     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
C     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
C     Y(I)=CZERO FOR I=NLAST+1,N
C
C***SEE ALSO  ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZUCHK, ZUNIK, ZUOIK
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZUNI1
C     COMPLEX CFN,CONE,CRSC,CSCL,CSR,CSS,CWRK,CZERO,C1,C2,PHI,RZ,SUM,S1,
C    *S2,Y,Z,ZETA1,ZETA2
      DOUBLE PRECISION ALIM, APHI, ASCLE, BRY, CONER, CRSC,
     * CSCL, CSRR, CSSR, CWRKI, CWRKR, C1R, C2I, C2M, C2R, ELIM, FN,
     * FNU, FNUL, PHII, PHIR, RAST, RS1, RZI, RZR, STI, STR, SUMI,
     * SUMR, S1I, S1R, S2I, S2R, TOL, YI, YR, ZEROI, ZEROR, ZETA1I,
     * ZETA1R, ZETA2I, ZETA2R, ZI, ZR, CYR, CYI, D1MACH, ZABS
      INTEGER I, IFLAG, INIT, K, KODE, M, N, ND, NLAST, NN, NUF, NW, NZ
      DIMENSION BRY(3), YR(N), YI(N), CWRKR(16), CWRKI(16), CSSR(3),
     * CSRR(3), CYR(2), CYI(2)
      EXTERNAL ZABS
      DATA ZEROR,ZEROI,CONER / 0.0D0, 0.0D0, 1.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZUNI1
      NZ = 0
      ND = N
      NLAST = 0
C-----------------------------------------------------------------------
C     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
C     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
C     EXP(ALIM)=EXP(ELIM)*TOL
C-----------------------------------------------------------------------
      CSCL = 1.0D0/TOL
      CRSC = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CRSC
      CSRR(1) = CRSC
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = 1.0D+3*D1MACH(1)/TOL
C-----------------------------------------------------------------------
C     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
C-----------------------------------------------------------------------
      FN = MAX(FNU,1.0D0)
      INIT = 0
      CALL ZUNIK(ZR, ZI, FN, 1, 1, TOL, INIT, PHIR, PHII, ZETA1R,
     * ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
      IF (KODE.EQ.1) GO TO 10
      STR = ZR + ZETA2R
      STI = ZI + ZETA2I
      RAST = FN/ZABS(STR,STI)
      STR = STR*RAST*RAST
      STI = -STI*RAST*RAST
      S1R = -ZETA1R + STR
      S1I = -ZETA1I + STI
      GO TO 20
   10 CONTINUE
      S1R = -ZETA1R + ZETA2R
      S1I = -ZETA1I + ZETA2I
   20 CONTINUE
      RS1 = S1R
      IF (ABS(RS1).GT.ELIM) GO TO 130
   30 CONTINUE
      NN = MIN(2,ND)
      DO 80 I=1,NN
        FN = FNU + (ND-I)
        INIT = 0
        CALL ZUNIK(ZR, ZI, FN, 1, 0, TOL, INIT, PHIR, PHII, ZETA1R,
     *   ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
        IF (KODE.EQ.1) GO TO 40
        STR = ZR + ZETA2R
        STI = ZI + ZETA2I
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = -ZETA1R + STR
        S1I = -ZETA1I + STI + ZI
        GO TO 50
   40   CONTINUE
        S1R = -ZETA1R + ZETA2R
        S1I = -ZETA1I + ZETA2I
   50   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = S1R
        IF (ABS(RS1).GT.ELIM) GO TO 110
        IF (I.EQ.1) IFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 60
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = ZABS(PHIR,PHII)
        RS1 = RS1 + LOG(APHI)
        IF (ABS(RS1).GT.ELIM) GO TO 110
        IF (I.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 60
        IF (I.EQ.1) IFLAG = 3
   60   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 IF ABS(S1).LT.ASCLE
C-----------------------------------------------------------------------
        S2R = PHIR*SUMR - PHII*SUMI
        S2I = PHIR*SUMI + PHII*SUMR
        STR = EXP(S1R)*CSSR(IFLAG)
        S1R = STR*COS(S1I)
        S1I = STR*SIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S2R*S1I + S2I*S1R
        S2R = STR
        IF (IFLAG.NE.1) GO TO 70
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 110
   70   CONTINUE
        CYR(I) = S2R
        CYI(I) = S2I
        M = ND - I + 1
        YR(M) = S2R*CSRR(IFLAG)
        YI(M) = S2I*CSRR(IFLAG)
   80 CONTINUE
      IF (ND.LE.2) GO TO 100
      RAST = 1.0D0/ZABS(ZR,ZI)
      STR = ZR*RAST
      STI = -ZI*RAST
      RZR = (STR+STR)*RAST
      RZI = (STI+STI)*RAST
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = D1MACH(2)
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = CSRR(IFLAG)
      ASCLE = BRY(IFLAG)
      K = ND - 2
      FN = K
      DO 90 I=3,ND
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FNU+FN)*(RZR*C2R-RZI*C2I)
        S2I = S1I + (FNU+FN)*(RZR*C2I+RZI*C2R)
        S1R = C2R
        S1I = C2I
        C2R = S2R*C1R
        C2I = S2I*C1R
        YR(K) = C2R
        YI(K) = C2I
        K = K - 1
        FN = FN - 1.0D0
        IF (IFLAG.GE.3) GO TO 90
        STR = ABS(C2R)
        STI = ABS(C2I)
        C2M = MAX(STR,STI)
        IF (C2M.LE.ASCLE) GO TO 90
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1R = S1R*C1R
        S1I = S1I*C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R*CSSR(IFLAG)
        S1I = S1I*CSSR(IFLAG)
        S2R = S2R*CSSR(IFLAG)
        S2I = S2I*CSSR(IFLAG)
        C1R = CSRR(IFLAG)
   90 CONTINUE
  100 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     SET UNDERFLOW AND UPDATE PARAMETERS
C-----------------------------------------------------------------------
  110 CONTINUE
      IF (RS1.GT.0.0D0) GO TO 120
      YR(ND) = ZEROR
      YI(ND) = ZEROI
      NZ = NZ + 1
      ND = ND - 1
      IF (ND.EQ.0) GO TO 100
      CALL ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM)
      IF (NUF.LT.0) GO TO 120
      ND = ND - NUF
      NZ = NZ + NUF
      IF (ND.EQ.0) GO TO 100
      FN = FNU + (ND-1)
      IF (FN.GE.FNUL) GO TO 30
      NLAST = ND
      RETURN
  120 CONTINUE
      NZ = -1
      RETURN
  130 CONTINUE
      IF (RS1.GT.0.0D0) GO TO 120
      NZ = N
      DO 140 I=1,N
        YR(I) = ZEROR
        YI(I) = ZEROI
  140 CONTINUE
      RETURN
      END
*DECK ZUNI2
      SUBROUTINE ZUNI2 (ZR, ZI, FNU, KODE, N, YR, YI, NZ, NLAST, FNUL,
     +   TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  ZUNI2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CUNI2-A, ZUNI2-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF
C     UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I
C     OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO.
C
C     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
C     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
C     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
C     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
C     Y(I)=CZERO FOR I=NLAST+1,N
C
C***SEE ALSO  ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZAIRY, ZUCHK, ZUNHJ, ZUOIK
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZUNI2
C     COMPLEX AI,ARG,ASUM,BSUM,CFN,CI,CID,CIP,CONE,CRSC,CSCL,CSR,CSS,
C    *CZERO,C1,C2,DAI,PHI,RZ,S1,S2,Y,Z,ZB,ZETA1,ZETA2,ZN
      DOUBLE PRECISION AARG, AIC, AII, AIR, ALIM, ANG, APHI, ARGI,
     * ARGR, ASCLE, ASUMI, ASUMR, BRY, BSUMI, BSUMR, CIDI, CIPI, CIPR,
     * CONER, CRSC, CSCL, CSRR, CSSR, C1R, C2I, C2M, C2R, DAII,
     * DAIR, ELIM, FN, FNU, FNUL, HPI, PHII, PHIR, RAST, RAZ, RS1, RZI,
     * RZR, STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR, ZBI, ZBR, ZEROI,
     * ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZI, ZNI, ZNR, ZR, CYR,
     * CYI, D1MACH, ZABS, CAR, SAR
      INTEGER I, IFLAG, IN, INU, J, K, KODE, N, NAI, ND, NDAI, NLAST,
     * NN, NUF, NW, NZ, IDUM
      DIMENSION BRY(3), YR(N), YI(N), CIPR(4), CIPI(4), CSSR(3),
     * CSRR(3), CYR(2), CYI(2)
      EXTERNAL ZABS
      DATA ZEROR,ZEROI,CONER / 0.0D0, 0.0D0, 1.0D0 /
      DATA CIPR(1),CIPI(1),CIPR(2),CIPI(2),CIPR(3),CIPI(3),CIPR(4),
     * CIPI(4)/ 1.0D0,0.0D0, 0.0D0,1.0D0, -1.0D0,0.0D0, 0.0D0,-1.0D0/
      DATA HPI, AIC  /
     1      1.57079632679489662D+00,     1.265512123484645396D+00/
C***FIRST EXECUTABLE STATEMENT  ZUNI2
      NZ = 0
      ND = N
      NLAST = 0
C-----------------------------------------------------------------------
C     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
C     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
C     EXP(ALIM)=EXP(ELIM)*TOL
C-----------------------------------------------------------------------
      CSCL = 1.0D0/TOL
      CRSC = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CRSC
      CSRR(1) = CRSC
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = 1.0D+3*D1MACH(1)/TOL
C-----------------------------------------------------------------------
C     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
C-----------------------------------------------------------------------
      ZNR = ZI
      ZNI = -ZR
      ZBR = ZR
      ZBI = ZI
      CIDI = -CONER
      INU = FNU
      ANG = HPI*(FNU-INU)
      C2R = COS(ANG)
      C2I = SIN(ANG)
      CAR = C2R
      SAR = C2I
      IN = INU + N - 1
      IN = MOD(IN,4) + 1
      STR = C2R*CIPR(IN) - C2I*CIPI(IN)
      C2I = C2R*CIPI(IN) + C2I*CIPR(IN)
      C2R = STR
      IF (ZI.GT.0.0D0) GO TO 10
      ZNR = -ZNR
      ZBI = -ZBI
      CIDI = -CIDI
      C2I = -C2I
   10 CONTINUE
C-----------------------------------------------------------------------
C     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
C-----------------------------------------------------------------------
      FN = MAX(FNU,1.0D0)
      CALL ZUNHJ(ZNR, ZNI, FN, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R,
     * ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
      IF (KODE.EQ.1) GO TO 20
      STR = ZBR + ZETA2R
      STI = ZBI + ZETA2I
      RAST = FN/ZABS(STR,STI)
      STR = STR*RAST*RAST
      STI = -STI*RAST*RAST
      S1R = -ZETA1R + STR
      S1I = -ZETA1I + STI
      GO TO 30
   20 CONTINUE
      S1R = -ZETA1R + ZETA2R
      S1I = -ZETA1I + ZETA2I
   30 CONTINUE
      RS1 = S1R
      IF (ABS(RS1).GT.ELIM) GO TO 150
   40 CONTINUE
      NN = MIN(2,ND)
      DO 90 I=1,NN
        FN = FNU + (ND-I)
        CALL ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIR, PHII, ARGR, ARGI,
     *   ZETA1R, ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
        IF (KODE.EQ.1) GO TO 50
        STR = ZBR + ZETA2R
        STI = ZBI + ZETA2I
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = -ZETA1R + STR
        S1I = -ZETA1I + STI + ABS(ZI)
        GO TO 60
   50   CONTINUE
        S1R = -ZETA1R + ZETA2R
        S1I = -ZETA1I + ZETA2I
   60   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = S1R
        IF (ABS(RS1).GT.ELIM) GO TO 120
        IF (I.EQ.1) IFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 70
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        APHI = ZABS(PHIR,PHII)
        AARG = ZABS(ARGR,ARGI)
        RS1 = RS1 + LOG(APHI) - 0.25D0*LOG(AARG) - AIC
        IF (ABS(RS1).GT.ELIM) GO TO 120
        IF (I.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 70
        IF (I.EQ.1) IFLAG = 3
   70   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
        CALL ZAIRY(ARGR, ARGI, 0, 2, AIR, AII, NAI, IDUM)
        CALL ZAIRY(ARGR, ARGI, 1, 2, DAIR, DAII, NDAI, IDUM)
        STR = DAIR*BSUMR - DAII*BSUMI
        STI = DAIR*BSUMI + DAII*BSUMR
        STR = STR + (AIR*ASUMR-AII*ASUMI)
        STI = STI + (AIR*ASUMI+AII*ASUMR)
        S2R = PHIR*STR - PHII*STI
        S2I = PHIR*STI + PHII*STR
        STR = EXP(S1R)*CSSR(IFLAG)
        S1R = STR*COS(S1I)
        S1I = STR*SIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S2R*S1I + S2I*S1R
        S2R = STR
        IF (IFLAG.NE.1) GO TO 80
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 120
   80   CONTINUE
        IF (ZI.LE.0.0D0) S2I = -S2I
        STR = S2R*C2R - S2I*C2I
        S2I = S2R*C2I + S2I*C2R
        S2R = STR
        CYR(I) = S2R
        CYI(I) = S2I
        J = ND - I + 1
        YR(J) = S2R*CSRR(IFLAG)
        YI(J) = S2I*CSRR(IFLAG)
        STR = -C2I*CIDI
        C2I = C2R*CIDI
        C2R = STR
   90 CONTINUE
      IF (ND.LE.2) GO TO 110
      RAZ = 1.0D0/ZABS(ZR,ZI)
      STR = ZR*RAZ
      STI = -ZI*RAZ
      RZR = (STR+STR)*RAZ
      RZI = (STI+STI)*RAZ
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = D1MACH(2)
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = CSRR(IFLAG)
      ASCLE = BRY(IFLAG)
      K = ND - 2
      FN = K
      DO 100 I=3,ND
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FNU+FN)*(RZR*C2R-RZI*C2I)
        S2I = S1I + (FNU+FN)*(RZR*C2I+RZI*C2R)
        S1R = C2R
        S1I = C2I
        C2R = S2R*C1R
        C2I = S2I*C1R
        YR(K) = C2R
        YI(K) = C2I
        K = K - 1
        FN = FN - 1.0D0
        IF (IFLAG.GE.3) GO TO 100
        STR = ABS(C2R)
        STI = ABS(C2I)
        C2M = MAX(STR,STI)
        IF (C2M.LE.ASCLE) GO TO 100
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1R = S1R*C1R
        S1I = S1I*C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R*CSSR(IFLAG)
        S1I = S1I*CSSR(IFLAG)
        S2R = S2R*CSSR(IFLAG)
        S2I = S2I*CSSR(IFLAG)
        C1R = CSRR(IFLAG)
  100 CONTINUE
  110 CONTINUE
      RETURN
  120 CONTINUE
      IF (RS1.GT.0.0D0) GO TO 140
C-----------------------------------------------------------------------
C     SET UNDERFLOW AND UPDATE PARAMETERS
C-----------------------------------------------------------------------
      YR(ND) = ZEROR
      YI(ND) = ZEROI
      NZ = NZ + 1
      ND = ND - 1
      IF (ND.EQ.0) GO TO 110
      CALL ZUOIK(ZR, ZI, FNU, KODE, 1, ND, YR, YI, NUF, TOL, ELIM, ALIM)
      IF (NUF.LT.0) GO TO 140
      ND = ND - NUF
      NZ = NZ + NUF
      IF (ND.EQ.0) GO TO 110
      FN = FNU + (ND-1)
      IF (FN.LT.FNUL) GO TO 130
C      FN = CIDI
C      J = NUF + 1
C      K = MOD(J,4) + 1
C      S1R = CIPR(K)
C      S1I = CIPI(K)
C      IF (FN.LT.0.0D0) S1I = -S1I
C      STR = C2R*S1R - C2I*S1I
C      C2I = C2R*S1I + C2I*S1R
C      C2R = STR
      IN = INU + ND - 1
      IN = MOD(IN,4) + 1
      C2R = CAR*CIPR(IN) - SAR*CIPI(IN)
      C2I = CAR*CIPI(IN) + SAR*CIPR(IN)
      IF (ZI.LE.0.0D0) C2I = -C2I
      GO TO 40
  130 CONTINUE
      NLAST = ND
      RETURN
  140 CONTINUE
      NZ = -1
      RETURN
  150 CONTINUE
      IF (RS1.GT.0.0D0) GO TO 140
      NZ = N
      DO 160 I=1,N
        YR(I) = ZEROR
        YI(I) = ZEROI
  160 CONTINUE
      RETURN
      END
*DECK ZUNIK
      SUBROUTINE ZUNIK (ZRR, ZRI, FNU, IKFLG, IPMTR, TOL, INIT, PHIR,
     +   PHII, ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
C***BEGIN PROLOGUE  ZUNIK
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CUNIK-A, ZUNIK-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C        ZUNIK COMPUTES PARAMETERS FOR THE UNIFORM ASYMPTOTIC
C        EXPANSIONS OF THE I AND K FUNCTIONS ON IKFLG= 1 OR 2
C        RESPECTIVELY BY
C
C        W(FNU,ZR) = PHI*EXP(ZETA)*SUM
C
C        WHERE       ZETA=-ZETA1 + ZETA2       OR
C                          ZETA1 - ZETA2
C
C        THE FIRST CALL MUST HAVE INIT=0. SUBSEQUENT CALLS WITH THE
C        SAME ZR AND FNU WILL RETURN THE I OR K FUNCTION ON IKFLG=
C        1 OR 2 WITH NO CHANGE IN INIT. CWRK IS A COMPLEX WORK
C        ARRAY. IPMTR=0 COMPUTES ALL PARAMETERS. IPMTR=1 COMPUTES PHI,
C        ZETA1,ZETA2.
C
C***SEE ALSO  ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, ZDIV, ZLOG, ZSQRT
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   930122  Added EXTERNAL statement with ZLOG and ZSQRT.  (RWC)
C***END PROLOGUE  ZUNIK
C     COMPLEX CFN,CON,CONE,CRFN,CWRK,CZERO,PHI,S,SR,SUM,T,T2,ZETA1,
C    *ZETA2,ZN,ZR
      DOUBLE PRECISION AC, C, CON, CONEI, CONER, CRFNI, CRFNR, CWRKI,
     * CWRKR, FNU, PHII, PHIR, RFN, SI, SR, SRI, SRR, STI, STR, SUMI,
     * SUMR, TEST, TI, TOL, TR, T2I, T2R, ZEROI, ZEROR, ZETA1I, ZETA1R,
     * ZETA2I, ZETA2R, ZNI, ZNR, ZRI, ZRR, D1MACH
      INTEGER I, IDUM, IKFLG, INIT, IPMTR, J, K, L
      DIMENSION C(120), CWRKR(16), CWRKI(16), CON(2)
      EXTERNAL ZLOG, ZSQRT
      DATA ZEROR,ZEROI,CONER,CONEI / 0.0D0, 0.0D0, 1.0D0, 0.0D0 /
      DATA CON(1), CON(2)  /
     1 3.98942280401432678D-01,  1.25331413731550025D+00 /
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2     C(19), C(20), C(21), C(22), C(23), C(24)/
     3     1.00000000000000000D+00,    -2.08333333333333333D-01,
     4     1.25000000000000000D-01,     3.34201388888888889D-01,
     5    -4.01041666666666667D-01,     7.03125000000000000D-02,
     6    -1.02581259645061728D+00,     1.84646267361111111D+00,
     7    -8.91210937500000000D-01,     7.32421875000000000D-02,
     8     4.66958442342624743D+00,    -1.12070026162229938D+01,
     9     8.78912353515625000D+00,    -2.36408691406250000D+00,
     A     1.12152099609375000D-01,    -2.82120725582002449D+01,
     B     8.46362176746007346D+01,    -9.18182415432400174D+01,
     C     4.25349987453884549D+01,    -7.36879435947963170D+00,
     D     2.27108001708984375D-01,     2.12570130039217123D+02,
     E    -7.65252468141181642D+02,     1.05999045252799988D+03/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3    -6.99579627376132541D+02,     2.18190511744211590D+02,
     4    -2.64914304869515555D+01,     5.72501420974731445D-01,
     5    -1.91945766231840700D+03,     8.06172218173730938D+03,
     6    -1.35865500064341374D+04,     1.16553933368645332D+04,
     7    -5.30564697861340311D+03,     1.20090291321635246D+03,
     8    -1.08090919788394656D+02,     1.72772750258445740D+00,
     9     2.02042913309661486D+04,    -9.69805983886375135D+04,
     A     1.92547001232531532D+05,    -2.03400177280415534D+05,
     B     1.22200464983017460D+05,    -4.11926549688975513D+04,
     C     7.10951430248936372D+03,    -4.93915304773088012D+02,
     D     6.07404200127348304D+00,    -2.42919187900551333D+05,
     E     1.31176361466297720D+06,    -2.99801591853810675D+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2     C(65), C(66), C(67), C(68), C(69), C(70), C(71), C(72)/
     3     3.76327129765640400D+06,    -2.81356322658653411D+06,
     4     1.26836527332162478D+06,    -3.31645172484563578D+05,
     5     4.52187689813627263D+04,    -2.49983048181120962D+03,
     6     2.43805296995560639D+01,     3.28446985307203782D+06,
     7    -1.97068191184322269D+07,     5.09526024926646422D+07,
     8    -7.41051482115326577D+07,     6.63445122747290267D+07,
     9    -3.75671766607633513D+07,     1.32887671664218183D+07,
     A    -2.78561812808645469D+06,     3.08186404612662398D+05,
     B    -1.38860897537170405D+04,     1.10017140269246738D+02,
     C    -4.93292536645099620D+07,     3.25573074185765749D+08,
     D    -9.39462359681578403D+08,     1.55359689957058006D+09,
     E    -1.62108055210833708D+09,     1.10684281682301447D+09/
      DATA C(73), C(74), C(75), C(76), C(77), C(78), C(79), C(80),
     1     C(81), C(82), C(83), C(84), C(85), C(86), C(87), C(88),
     2     C(89), C(90), C(91), C(92), C(93), C(94), C(95), C(96)/
     3    -4.95889784275030309D+08,     1.42062907797533095D+08,
     4    -2.44740627257387285D+07,     2.24376817792244943D+06,
     5    -8.40054336030240853D+04,     5.51335896122020586D+02,
     6     8.14789096118312115D+08,    -5.86648149205184723D+09,
     7     1.86882075092958249D+10,    -3.46320433881587779D+10,
     8     4.12801855797539740D+10,    -3.30265997498007231D+10,
     9     1.79542137311556001D+10,    -6.56329379261928433D+09,
     A     1.55927986487925751D+09,    -2.25105661889415278D+08,
     B     1.73951075539781645D+07,    -5.49842327572288687D+05,
     C     3.03809051092238427D+03,    -1.46792612476956167D+10,
     D     1.14498237732025810D+11,    -3.99096175224466498D+11,
     E     8.19218669548577329D+11,    -1.09837515608122331D+12/
      DATA C(97), C(98), C(99), C(100), C(101), C(102), C(103), C(104),
     1     C(105), C(106), C(107), C(108), C(109), C(110), C(111),
     2     C(112), C(113), C(114), C(115), C(116), C(117), C(118)/
     3     1.00815810686538209D+12,    -6.45364869245376503D+11,
     4     2.87900649906150589D+11,    -8.78670721780232657D+10,
     5     1.76347306068349694D+10,    -2.16716498322379509D+09,
     6     1.43157876718888981D+08,    -3.87183344257261262D+06,
     7     1.82577554742931747D+04,     2.86464035717679043D+11,
     8    -2.40629790002850396D+12,     9.10934118523989896D+12,
     9    -2.05168994109344374D+13,     3.05651255199353206D+13,
     A    -3.16670885847851584D+13,     2.33483640445818409D+13,
     B    -1.23204913055982872D+13,     4.61272578084913197D+12,
     C    -1.19655288019618160D+12,     2.05914503232410016D+11,
     D    -2.18229277575292237D+10,     1.24700929351271032D+09/
      DATA C(119), C(120)/
     1    -2.91883881222208134D+07,     1.18838426256783253D+05/
C***FIRST EXECUTABLE STATEMENT  ZUNIK
      IF (INIT.NE.0) GO TO 40
C-----------------------------------------------------------------------
C     INITIALIZE ALL VARIABLES
C-----------------------------------------------------------------------
      RFN = 1.0D0/FNU
C-----------------------------------------------------------------------
C     OVERFLOW TEST (ZR/FNU TOO SMALL)
C-----------------------------------------------------------------------
      TEST = D1MACH(1)*1.0D+3
      AC = FNU*TEST
      IF (ABS(ZRR).GT.AC .OR. ABS(ZRI).GT.AC) GO TO 15
      ZETA1R = 2.0D0*ABS(LOG(TEST))+FNU
      ZETA1I = 0.0D0
      ZETA2R = FNU
      ZETA2I = 0.0D0
      PHIR = 1.0D0
      PHII = 0.0D0
      RETURN
   15 CONTINUE
      TR = ZRR*RFN
      TI = ZRI*RFN
      SR = CONER + (TR*TR-TI*TI)
      SI = CONEI + (TR*TI+TI*TR)
      CALL ZSQRT(SR, SI, SRR, SRI)
      STR = CONER + SRR
      STI = CONEI + SRI
      CALL ZDIV(STR, STI, TR, TI, ZNR, ZNI)
      CALL ZLOG(ZNR, ZNI, STR, STI, IDUM)
      ZETA1R = FNU*STR
      ZETA1I = FNU*STI
      ZETA2R = FNU*SRR
      ZETA2I = FNU*SRI
      CALL ZDIV(CONER, CONEI, SRR, SRI, TR, TI)
      SRR = TR*RFN
      SRI = TI*RFN
      CALL ZSQRT(SRR, SRI, CWRKR(16), CWRKI(16))
      PHIR = CWRKR(16)*CON(IKFLG)
      PHII = CWRKI(16)*CON(IKFLG)
      IF (IPMTR.NE.0) RETURN
      CALL ZDIV(CONER, CONEI, SR, SI, T2R, T2I)
      CWRKR(1) = CONER
      CWRKI(1) = CONEI
      CRFNR = CONER
      CRFNI = CONEI
      AC = 1.0D0
      L = 1
      DO 20 K=2,15
        SR = ZEROR
        SI = ZEROI
        DO 10 J=1,K
          L = L + 1
          STR = SR*T2R - SI*T2I + C(L)
          SI = SR*T2I + SI*T2R
          SR = STR
   10   CONTINUE
        STR = CRFNR*SRR - CRFNI*SRI
        CRFNI = CRFNR*SRI + CRFNI*SRR
        CRFNR = STR
        CWRKR(K) = CRFNR*SR - CRFNI*SI
        CWRKI(K) = CRFNR*SI + CRFNI*SR
        AC = AC*RFN
        TEST = ABS(CWRKR(K)) + ABS(CWRKI(K))
        IF (AC.LT.TOL .AND. TEST.LT.TOL) GO TO 30
   20 CONTINUE
      K = 15
   30 CONTINUE
      INIT = K
   40 CONTINUE
      IF (IKFLG.EQ.2) GO TO 60
C-----------------------------------------------------------------------
C     COMPUTE SUM FOR THE I FUNCTION
C-----------------------------------------------------------------------
      SR = ZEROR
      SI = ZEROI
      DO 50 I=1,INIT
        SR = SR + CWRKR(I)
        SI = SI + CWRKI(I)
   50 CONTINUE
      SUMR = SR
      SUMI = SI
      PHIR = CWRKR(16)*CON(1)
      PHII = CWRKI(16)*CON(1)
      RETURN
   60 CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE SUM FOR THE K FUNCTION
C-----------------------------------------------------------------------
      SR = ZEROR
      SI = ZEROI
      TR = CONER
      DO 70 I=1,INIT
        SR = SR + TR*CWRKR(I)
        SI = SI + TR*CWRKI(I)
        TR = -TR
   70 CONTINUE
      SUMR = SR
      SUMI = SI
      PHIR = CWRKR(16)*CON(2)
      PHII = CWRKI(16)*CON(2)
      RETURN
      END
*DECK ZUNK1
      SUBROUTINE ZUNK1 (ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
     +   ALIM)
C***BEGIN PROLOGUE  ZUNK1
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CUNK1-A, ZUNK1-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
C     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
C     UNIFORM ASYMPTOTIC EXPANSION.
C     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
C     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
C
C***SEE ALSO  ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZS1S2, ZUCHK, ZUNIK
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZUNK1
C     COMPLEX CFN,CK,CONE,CRSC,CS,CSCL,CSGN,CSPN,CSR,CSS,CWRK,CY,CZERO,
C    *C1,C2,PHI,PHID,RZ,SUM,SUMD,S1,S2,Y,Z,ZETA1,ZETA1D,ZETA2,ZETA2D,ZR
      DOUBLE PRECISION ALIM, ANG, APHI, ASC, ASCLE, BRY, CKI, CKR,
     * CONER, CRSC, CSCL, CSGNI, CSPNI, CSPNR, CSR, CSRR, CSSR,
     * CWRKI, CWRKR, CYI, CYR, C1I, C1R, C2I, C2M, C2R, ELIM, FMR, FN,
     * FNF, FNU, PHIDI, PHIDR, PHII, PHIR, PI, RAST, RAZR, RS1, RZI,
     * RZR, SGN, STI, STR, SUMDI, SUMDR, SUMI, SUMR, S1I, S1R, S2I,
     * S2R, TOL, YI, YR, ZEROI, ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R,
     * ZET1DI, ZET1DR, ZET2DI, ZET2DR, ZI, ZR, ZRI, ZRR, D1MACH, ZABS
      INTEGER I, IB, IFLAG, IFN, IL, INIT, INU, IUF, K, KDFLG, KFLAG,
     * KK, KODE, MR, N, NW, NZ, INITD, IC, IPARD, J, M
      DIMENSION BRY(3), INIT(2), YR(N), YI(N), SUMR(2), SUMI(2),
     * ZETA1R(2), ZETA1I(2), ZETA2R(2), ZETA2I(2), CYR(2), CYI(2),
     * CWRKR(16,3), CWRKI(16,3), CSSR(3), CSRR(3), PHIR(2), PHII(2)
      EXTERNAL ZABS
      DATA ZEROR,ZEROI,CONER / 0.0D0, 0.0D0, 1.0D0 /
      DATA PI / 3.14159265358979324D0 /
C***FIRST EXECUTABLE STATEMENT  ZUNK1
      KDFLG = 1
      NZ = 0
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
C     THE UNDERFLOW LIMIT
C-----------------------------------------------------------------------
      CSCL = 1.0D0/TOL
      CRSC = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CRSC
      CSRR(1) = CRSC
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = 1.0D+3*D1MACH(1)/TOL
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = D1MACH(2)
      ZRR = ZR
      ZRI = ZI
      IF (ZR.GE.0.0D0) GO TO 10
      ZRR = -ZR
      ZRI = -ZI
   10 CONTINUE
      J = 2
      DO 70 I=1,N
C-----------------------------------------------------------------------
C     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
C-----------------------------------------------------------------------
        J = 3 - J
        FN = FNU + (I-1)
        INIT(J) = 0
        CALL ZUNIK(ZRR, ZRI, FN, 2, 0, TOL, INIT(J), PHIR(J), PHII(J),
     *   ZETA1R(J), ZETA1I(J), ZETA2R(J), ZETA2I(J), SUMR(J), SUMI(J),
     *   CWRKR(1,J), CWRKI(1,J))
        IF (KODE.EQ.1) GO TO 20
        STR = ZRR + ZETA2R(J)
        STI = ZRI + ZETA2I(J)
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = ZETA1R(J) - STR
        S1I = ZETA1I(J) - STI
        GO TO 30
   20   CONTINUE
        S1R = ZETA1R(J) - ZETA2R(J)
        S1I = ZETA1I(J) - ZETA2I(J)
   30   CONTINUE
        RS1 = S1R
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        IF (ABS(RS1).GT.ELIM) GO TO 60
        IF (KDFLG.EQ.1) KFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 40
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = ZABS(PHIR(J),PHII(J))
        RS1 = RS1 + LOG(APHI)
        IF (ABS(RS1).GT.ELIM) GO TO 60
        IF (KDFLG.EQ.1) KFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 40
        IF (KDFLG.EQ.1) KFLAG = 3
   40   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
        S2R = PHIR(J)*SUMR(J) - PHII(J)*SUMI(J)
        S2I = PHIR(J)*SUMI(J) + PHII(J)*SUMR(J)
        STR = EXP(S1R)*CSSR(KFLAG)
        S1R = STR*COS(S1I)
        S1I = STR*SIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S1R*S2I + S2R*S1I
        S2R = STR
        IF (KFLAG.NE.1) GO TO 50
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 60
   50   CONTINUE
        CYR(KDFLG) = S2R
        CYI(KDFLG) = S2I
        YR(I) = S2R*CSRR(KFLAG)
        YI(I) = S2I*CSRR(KFLAG)
        IF (KDFLG.EQ.2) GO TO 75
        KDFLG = 2
        GO TO 70
   60   CONTINUE
        IF (RS1.GT.0.0D0) GO TO 300
C-----------------------------------------------------------------------
C     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
        IF (ZR.LT.0.0D0) GO TO 300
        KDFLG = 1
        YR(I)=ZEROR
        YI(I)=ZEROI
        NZ=NZ+1
        IF (I.EQ.1) GO TO 70
        IF ((YR(I-1).EQ.ZEROR).AND.(YI(I-1).EQ.ZEROI)) GO TO 70
        YR(I-1)=ZEROR
        YI(I-1)=ZEROI
        NZ=NZ+1
   70 CONTINUE
      I = N
   75 CONTINUE
      RAZR = 1.0D0/ZABS(ZRR,ZRI)
      STR = ZRR*RAZR
      STI = -ZRI*RAZR
      RZR = (STR+STR)*RAZR
      RZI = (STI+STI)*RAZR
      CKR = FN*RZR
      CKI = FN*RZI
      IB = I + 1
      IF (N.LT.IB) GO TO 160
C-----------------------------------------------------------------------
C     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO
C     ON UNDERFLOW.
C-----------------------------------------------------------------------
      FN = FNU + (N-1)
      IPARD = 1
      IF (MR.NE.0) IPARD = 0
      INITD = 0
      CALL ZUNIK(ZRR, ZRI, FN, 2, IPARD, TOL, INITD, PHIDR, PHIDI,
     * ZET1DR, ZET1DI, ZET2DR, ZET2DI, SUMDR, SUMDI, CWRKR(1,3),
     * CWRKI(1,3))
      IF (KODE.EQ.1) GO TO 80
      STR = ZRR + ZET2DR
      STI = ZRI + ZET2DI
      RAST = FN/ZABS(STR,STI)
      STR = STR*RAST*RAST
      STI = -STI*RAST*RAST
      S1R = ZET1DR - STR
      S1I = ZET1DI - STI
      GO TO 90
   80 CONTINUE
      S1R = ZET1DR - ZET2DR
      S1I = ZET1DI - ZET2DI
   90 CONTINUE
      RS1 = S1R
      IF (ABS(RS1).GT.ELIM) GO TO 95
      IF (ABS(RS1).LT.ALIM) GO TO 100
C-----------------------------------------------------------------------
C     REFINE ESTIMATE AND TEST
C-----------------------------------------------------------------------
      APHI = ZABS(PHIDR,PHIDI)
      RS1 = RS1+LOG(APHI)
      IF (ABS(RS1).LT.ELIM) GO TO 100
   95 CONTINUE
      IF (ABS(RS1).GT.0.0D0) GO TO 300
C-----------------------------------------------------------------------
C     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
      IF (ZR.LT.0.0D0) GO TO 300
      NZ = N
      DO 96 I=1,N
        YR(I) = ZEROR
        YI(I) = ZEROI
   96 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     FORWARD RECUR FOR REMAINDER OF THE SEQUENCE
C-----------------------------------------------------------------------
  100 CONTINUE
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = CSRR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 120 I=IB,N
        C2R = S2R
        C2I = S2I
        S2R = CKR*C2R - CKI*C2I + S1R
        S2I = CKR*C2I + CKI*C2R + S1I
        S1R = C2R
        S1I = C2I
        CKR = CKR + RZR
        CKI = CKI + RZI
        C2R = S2R*C1R
        C2I = S2I*C1R
        YR(I) = C2R
        YI(I) = C2I
        IF (KFLAG.GE.3) GO TO 120
        STR = ABS(C2R)
        STI = ABS(C2I)
        C2M = MAX(STR,STI)
        IF (C2M.LE.ASCLE) GO TO 120
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1R = S1R*C1R
        S1I = S1I*C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R*CSSR(KFLAG)
        S1I = S1I*CSSR(KFLAG)
        S2R = S2R*CSSR(KFLAG)
        S2I = S2I*CSSR(KFLAG)
        C1R = CSRR(KFLAG)
  120 CONTINUE
  160 CONTINUE
      IF (MR.EQ.0) RETURN
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0D0
C-----------------------------------------------------------------------
      NZ = 0
      FMR = MR
      SGN = -DSIGN(PI,FMR)
C-----------------------------------------------------------------------
C     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
C-----------------------------------------------------------------------
      CSGNI = SGN
      INU = FNU
      FNF = FNU - INU
      IFN = INU + N - 1
      ANG = FNF*SGN
      CSPNR = COS(ANG)
      CSPNI = SIN(ANG)
      IF (MOD(IFN,2).EQ.0) GO TO 170
      CSPNR = -CSPNR
      CSPNI = -CSPNI
  170 CONTINUE
      ASC = BRY(1)
      IUF = 0
      KK = N
      KDFLG = 1
      IB = IB - 1
      IC = IB - 1
      DO 270 K=1,N
        FN = FNU + (KK-1)
C-----------------------------------------------------------------------
C     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
C     FUNCTION ABOVE
C-----------------------------------------------------------------------
        M=3
        IF (N.GT.2) GO TO 175
  172   CONTINUE
        INITD = INIT(J)
        PHIDR = PHIR(J)
        PHIDI = PHII(J)
        ZET1DR = ZETA1R(J)
        ZET1DI = ZETA1I(J)
        ZET2DR = ZETA2R(J)
        ZET2DI = ZETA2I(J)
        SUMDR = SUMR(J)
        SUMDI = SUMI(J)
        M = J
        J = 3 - J
        GO TO 180
  175   CONTINUE
        IF ((KK.EQ.N).AND.(IB.LT.N)) GO TO 180
        IF ((KK.EQ.IB).OR.(KK.EQ.IC)) GO TO 172
        INITD = 0
  180   CONTINUE
        CALL ZUNIK(ZRR, ZRI, FN, 1, 0, TOL, INITD, PHIDR, PHIDI,
     *   ZET1DR, ZET1DI, ZET2DR, ZET2DI, SUMDR, SUMDI,
     *   CWRKR(1,M), CWRKI(1,M))
        IF (KODE.EQ.1) GO TO 200
        STR = ZRR + ZET2DR
        STI = ZRI + ZET2DI
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = -ZET1DR + STR
        S1I = -ZET1DI + STI
        GO TO 210
  200   CONTINUE
        S1R = -ZET1DR + ZET2DR
        S1I = -ZET1DI + ZET2DI
  210   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = S1R
        IF (ABS(RS1).GT.ELIM) GO TO 260
        IF (KDFLG.EQ.1) IFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 220
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = ZABS(PHIDR,PHIDI)
        RS1 = RS1 + LOG(APHI)
        IF (ABS(RS1).GT.ELIM) GO TO 260
        IF (KDFLG.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 220
        IF (KDFLG.EQ.1) IFLAG = 3
  220   CONTINUE
        STR = PHIDR*SUMDR - PHIDI*SUMDI
        STI = PHIDR*SUMDI + PHIDI*SUMDR
        S2R = -CSGNI*STI
        S2I = CSGNI*STR
        STR = EXP(S1R)*CSSR(IFLAG)
        S1R = STR*COS(S1I)
        S1I = STR*SIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S2R*S1I + S2I*S1R
        S2R = STR
        IF (IFLAG.NE.1) GO TO 230
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.EQ.0) GO TO 230
        S2R = ZEROR
        S2I = ZEROI
  230   CONTINUE
        CYR(KDFLG) = S2R
        CYI(KDFLG) = S2I
        C2R = S2R
        C2I = S2I
        S2R = S2R*CSRR(IFLAG)
        S2I = S2I*CSRR(IFLAG)
C-----------------------------------------------------------------------
C     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
C-----------------------------------------------------------------------
        S1R = YR(KK)
        S1I = YI(KK)
        IF (KODE.EQ.1) GO TO 250
        CALL ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  250   CONTINUE
        YR(KK) = S1R*CSPNR - S1I*CSPNI + S2R
        YI(KK) = CSPNR*S1I + CSPNI*S1R + S2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        IF (C2R.NE.0.0D0 .OR. C2I.NE.0.0D0) GO TO 255
        KDFLG = 1
        GO TO 270
  255   CONTINUE
        IF (KDFLG.EQ.2) GO TO 275
        KDFLG = 2
        GO TO 270
  260   CONTINUE
        IF (RS1.GT.0.0D0) GO TO 300
        S2R = ZEROR
        S2I = ZEROI
        GO TO 230
  270 CONTINUE
      K = N
  275 CONTINUE
      IL = N - K
      IF (IL.EQ.0) RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
C     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
C     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
C-----------------------------------------------------------------------
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      CSR = CSRR(IFLAG)
      ASCLE = BRY(IFLAG)
      FN = INU+IL
      DO 290 I=1,IL
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FN+FNF)*(RZR*C2R-RZI*C2I)
        S2I = S1I + (FN+FNF)*(RZR*C2I+RZI*C2R)
        S1R = C2R
        S1I = C2I
        FN = FN - 1.0D0
        C2R = S2R*CSR
        C2I = S2I*CSR
        CKR = C2R
        CKI = C2I
        C1R = YR(KK)
        C1I = YI(KK)
        IF (KODE.EQ.1) GO TO 280
        CALL ZS1S2(ZRR, ZRI, C1R, C1I, C2R, C2I, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  280   CONTINUE
        YR(KK) = C1R*CSPNR - C1I*CSPNI + C2R
        YI(KK) = C1R*CSPNI + C1I*CSPNR + C2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        IF (IFLAG.GE.3) GO TO 290
        C2R = ABS(CKR)
        C2I = ABS(CKI)
        C2M = MAX(C2R,C2I)
        IF (C2M.LE.ASCLE) GO TO 290
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1R = S1R*CSR
        S1I = S1I*CSR
        S2R = CKR
        S2I = CKI
        S1R = S1R*CSSR(IFLAG)
        S1I = S1I*CSSR(IFLAG)
        S2R = S2R*CSSR(IFLAG)
        S2I = S2I*CSSR(IFLAG)
        CSR = CSRR(IFLAG)
  290 CONTINUE
      RETURN
  300 CONTINUE
      NZ = -1
      RETURN
      END
*DECK ZUNK2
      SUBROUTINE ZUNK2 (ZR, ZI, FNU, KODE, MR, N, YR, YI, NZ, TOL, ELIM,
     +   ALIM)
C***BEGIN PROLOGUE  ZUNK2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CUNK2-A, ZUNK2-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZUNK2 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
C     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
C     UNIFORM ASYMPTOTIC EXPANSIONS FOR H(KIND,FNU,ZN) AND J(FNU,ZN)
C     WHERE ZN IS IN THE RIGHT HALF PLANE, KIND=(3-MR)/2, MR=+1 OR
C     -1. HERE ZN=ZR*I OR -ZR*I WHERE ZR=Z IF Z IS IN THE RIGHT
C     HALF PLANE OR ZR=-Z IF Z IS IN THE LEFT HALF PLANE. MR INDIC-
C     ATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
C     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
C
C***SEE ALSO  ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZAIRY, ZS1S2, ZUCHK, ZUNHJ
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZUNK2
C     COMPLEX AI,ARG,ARGD,ASUM,ASUMD,BSUM,BSUMD,CFN,CI,CIP,CK,CONE,CRSC,
C    *CR1,CR2,CS,CSCL,CSGN,CSPN,CSR,CSS,CY,CZERO,C1,C2,DAI,PHI,PHID,RZ,
C    *S1,S2,Y,Z,ZB,ZETA1,ZETA1D,ZETA2,ZETA2D,ZN,ZR
      DOUBLE PRECISION AARG, AIC, AII, AIR, ALIM, ANG, APHI, ARGDI,
     * ARGDR, ARGI, ARGR, ASC, ASCLE, ASUMDI, ASUMDR, ASUMI, ASUMR,
     * BRY, BSUMDI, BSUMDR, BSUMI, BSUMR, CAR, CIPI, CIPR, CKI, CKR,
     * CONER, CRSC, CR1I, CR1R, CR2I, CR2R, CSCL, CSGNI, CSI,
     * CSPNI, CSPNR, CSR, CSRR, CSSR, CYI, CYR, C1I, C1R, C2I, C2M,
     * C2R, DAII, DAIR, ELIM, FMR, FN, FNF, FNU, HPI, PHIDI, PHIDR,
     * PHII, PHIR, PI, PTI, PTR, RAST, RAZR, RS1, RZI, RZR, SAR, SGN,
     * STI, STR, S1I, S1R, S2I, S2R, TOL, YI, YR, YY, ZBI, ZBR, ZEROI,
     * ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZET1DI, ZET1DR, ZET2DI,
     * ZET2DR, ZI, ZNI, ZNR, ZR, ZRI, ZRR, D1MACH, ZABS
      INTEGER I, IB, IFLAG, IFN, IL, IN, INU, IUF, K, KDFLG, KFLAG, KK,
     * KODE, MR, N, NAI, NDAI, NW, NZ, IDUM, J, IPARD, IC
      DIMENSION BRY(3), YR(N), YI(N), ASUMR(2), ASUMI(2), BSUMR(2),
     * BSUMI(2), PHIR(2), PHII(2), ARGR(2), ARGI(2), ZETA1R(2),
     * ZETA1I(2), ZETA2R(2), ZETA2I(2), CYR(2), CYI(2), CIPR(4),
     * CIPI(4), CSSR(3), CSRR(3)
      EXTERNAL ZABS
      DATA ZEROR,ZEROI,CONER,CR1R,CR1I,CR2R,CR2I /
     1         0.0D0, 0.0D0, 1.0D0,
     1 1.0D0,1.73205080756887729D0 , -0.5D0,-8.66025403784438647D-01 /
      DATA HPI, PI, AIC /
     1     1.57079632679489662D+00,     3.14159265358979324D+00,
     1     1.26551212348464539D+00/
      DATA CIPR(1),CIPI(1),CIPR(2),CIPI(2),CIPR(3),CIPI(3),CIPR(4),
     * CIPI(4) /
     1  1.0D0,0.0D0 ,  0.0D0,-1.0D0 ,  -1.0D0,0.0D0 ,  0.0D0,1.0D0 /
C***FIRST EXECUTABLE STATEMENT  ZUNK2
      KDFLG = 1
      NZ = 0
C-----------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
C     THE UNDERFLOW LIMIT
C-----------------------------------------------------------------------
      CSCL = 1.0D0/TOL
      CRSC = TOL
      CSSR(1) = CSCL
      CSSR(2) = CONER
      CSSR(3) = CRSC
      CSRR(1) = CRSC
      CSRR(2) = CONER
      CSRR(3) = CSCL
      BRY(1) = 1.0D+3*D1MACH(1)/TOL
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = D1MACH(2)
      ZRR = ZR
      ZRI = ZI
      IF (ZR.GE.0.0D0) GO TO 10
      ZRR = -ZR
      ZRI = -ZI
   10 CONTINUE
      YY = ZRI
      ZNR = ZRI
      ZNI = -ZRR
      ZBR = ZRR
      ZBI = ZRI
      INU = FNU
      FNF = FNU - INU
      ANG = -HPI*FNF
      CAR = COS(ANG)
      SAR = SIN(ANG)
      C2R = HPI*SAR
      C2I = -HPI*CAR
      KK = MOD(INU,4) + 1
      STR = C2R*CIPR(KK) - C2I*CIPI(KK)
      STI = C2R*CIPI(KK) + C2I*CIPR(KK)
      CSR = CR1R*STR - CR1I*STI
      CSI = CR1R*STI + CR1I*STR
      IF (YY.GT.0.0D0) GO TO 20
      ZNR = -ZNR
      ZBI = -ZBI
   20 CONTINUE
C-----------------------------------------------------------------------
C     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST
C     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY
C     CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS
C-----------------------------------------------------------------------
      J = 2
      DO 80 I=1,N
C-----------------------------------------------------------------------
C     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
C-----------------------------------------------------------------------
        J = 3 - J
        FN = FNU + (I-1)
        CALL ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIR(J), PHII(J), ARGR(J),
     *   ARGI(J), ZETA1R(J), ZETA1I(J), ZETA2R(J), ZETA2I(J), ASUMR(J),
     *   ASUMI(J), BSUMR(J), BSUMI(J))
        IF (KODE.EQ.1) GO TO 30
        STR = ZBR + ZETA2R(J)
        STI = ZBI + ZETA2I(J)
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = ZETA1R(J) - STR
        S1I = ZETA1I(J) - STI
        GO TO 40
   30   CONTINUE
        S1R = ZETA1R(J) - ZETA2R(J)
        S1I = ZETA1I(J) - ZETA2I(J)
   40   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = S1R
        IF (ABS(RS1).GT.ELIM) GO TO 70
        IF (KDFLG.EQ.1) KFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 50
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = ZABS(PHIR(J),PHII(J))
        AARG = ZABS(ARGR(J),ARGI(J))
        RS1 = RS1 + LOG(APHI) - 0.25D0*LOG(AARG) - AIC
        IF (ABS(RS1).GT.ELIM) GO TO 70
        IF (KDFLG.EQ.1) KFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 50
        IF (KDFLG.EQ.1) KFLAG = 3
   50   CONTINUE
C-----------------------------------------------------------------------
C     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C     EXPONENT EXTREMES
C-----------------------------------------------------------------------
        C2R = ARGR(J)*CR2R - ARGI(J)*CR2I
        C2I = ARGR(J)*CR2I + ARGI(J)*CR2R
        CALL ZAIRY(C2R, C2I, 0, 2, AIR, AII, NAI, IDUM)
        CALL ZAIRY(C2R, C2I, 1, 2, DAIR, DAII, NDAI, IDUM)
        STR = DAIR*BSUMR(J) - DAII*BSUMI(J)
        STI = DAIR*BSUMI(J) + DAII*BSUMR(J)
        PTR = STR*CR2R - STI*CR2I
        PTI = STR*CR2I + STI*CR2R
        STR = PTR + (AIR*ASUMR(J)-AII*ASUMI(J))
        STI = PTI + (AIR*ASUMI(J)+AII*ASUMR(J))
        PTR = STR*PHIR(J) - STI*PHII(J)
        PTI = STR*PHII(J) + STI*PHIR(J)
        S2R = PTR*CSR - PTI*CSI
        S2I = PTR*CSI + PTI*CSR
        STR = EXP(S1R)*CSSR(KFLAG)
        S1R = STR*COS(S1I)
        S1I = STR*SIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S1R*S2I + S2R*S1I
        S2R = STR
        IF (KFLAG.NE.1) GO TO 60
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.NE.0) GO TO 70
   60   CONTINUE
        IF (YY.LE.0.0D0) S2I = -S2I
        CYR(KDFLG) = S2R
        CYI(KDFLG) = S2I
        YR(I) = S2R*CSRR(KFLAG)
        YI(I) = S2I*CSRR(KFLAG)
        STR = CSI
        CSI = -CSR
        CSR = STR
        IF (KDFLG.EQ.2) GO TO 85
        KDFLG = 2
        GO TO 80
   70   CONTINUE
        IF (RS1.GT.0.0D0) GO TO 320
C-----------------------------------------------------------------------
C     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
        IF (ZR.LT.0.0D0) GO TO 320
        KDFLG = 1
        YR(I)=ZEROR
        YI(I)=ZEROI
        NZ=NZ+1
        STR = CSI
        CSI =-CSR
        CSR = STR
        IF (I.EQ.1) GO TO 80
        IF ((YR(I-1).EQ.ZEROR).AND.(YI(I-1).EQ.ZEROI)) GO TO 80
        YR(I-1)=ZEROR
        YI(I-1)=ZEROI
        NZ=NZ+1
   80 CONTINUE
      I = N
   85 CONTINUE
      RAZR = 1.0D0/ZABS(ZRR,ZRI)
      STR = ZRR*RAZR
      STI = -ZRI*RAZR
      RZR = (STR+STR)*RAZR
      RZI = (STI+STI)*RAZR
      CKR = FN*RZR
      CKI = FN*RZI
      IB = I + 1
      IF (N.LT.IB) GO TO 180
C-----------------------------------------------------------------------
C     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO
C     ON UNDERFLOW.
C-----------------------------------------------------------------------
      FN = FNU + (N-1)
      IPARD = 1
      IF (MR.NE.0) IPARD = 0
      CALL ZUNHJ(ZNR, ZNI, FN, IPARD, TOL, PHIDR, PHIDI, ARGDR, ARGDI,
     * ZET1DR, ZET1DI, ZET2DR, ZET2DI, ASUMDR, ASUMDI, BSUMDR, BSUMDI)
      IF (KODE.EQ.1) GO TO 90
      STR = ZBR + ZET2DR
      STI = ZBI + ZET2DI
      RAST = FN/ZABS(STR,STI)
      STR = STR*RAST*RAST
      STI = -STI*RAST*RAST
      S1R = ZET1DR - STR
      S1I = ZET1DI - STI
      GO TO 100
   90 CONTINUE
      S1R = ZET1DR - ZET2DR
      S1I = ZET1DI - ZET2DI
  100 CONTINUE
      RS1 = S1R
      IF (ABS(RS1).GT.ELIM) GO TO 105
      IF (ABS(RS1).LT.ALIM) GO TO 120
C-----------------------------------------------------------------------
C     REFINE ESTIMATE AND TEST
C-----------------------------------------------------------------------
      APHI = ZABS(PHIDR,PHIDI)
      RS1 = RS1+LOG(APHI)
      IF (ABS(RS1).LT.ELIM) GO TO 120
  105 CONTINUE
      IF (RS1.GT.0.0D0) GO TO 320
C-----------------------------------------------------------------------
C     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C-----------------------------------------------------------------------
      IF (ZR.LT.0.0D0) GO TO 320
      NZ = N
      DO 106 I=1,N
        YR(I) = ZEROR
        YI(I) = ZEROI
  106 CONTINUE
      RETURN
  120 CONTINUE
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      C1R = CSRR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 130 I=IB,N
        C2R = S2R
        C2I = S2I
        S2R = CKR*C2R - CKI*C2I + S1R
        S2I = CKR*C2I + CKI*C2R + S1I
        S1R = C2R
        S1I = C2I
        CKR = CKR + RZR
        CKI = CKI + RZI
        C2R = S2R*C1R
        C2I = S2I*C1R
        YR(I) = C2R
        YI(I) = C2I
        IF (KFLAG.GE.3) GO TO 130
        STR = ABS(C2R)
        STI = ABS(C2I)
        C2M = MAX(STR,STI)
        IF (C2M.LE.ASCLE) GO TO 130
        KFLAG = KFLAG + 1
        ASCLE = BRY(KFLAG)
        S1R = S1R*C1R
        S1I = S1I*C1R
        S2R = C2R
        S2I = C2I
        S1R = S1R*CSSR(KFLAG)
        S1I = S1I*CSSR(KFLAG)
        S2R = S2R*CSSR(KFLAG)
        S2I = S2I*CSSR(KFLAG)
        C1R = CSRR(KFLAG)
  130 CONTINUE
  180 CONTINUE
      IF (MR.EQ.0) RETURN
C-----------------------------------------------------------------------
C     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0D0
C-----------------------------------------------------------------------
      NZ = 0
      FMR = MR
      SGN = -DSIGN(PI,FMR)
C-----------------------------------------------------------------------
C     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
C-----------------------------------------------------------------------
      CSGNI = SGN
      IF (YY.LE.0.0D0) CSGNI = -CSGNI
      IFN = INU + N - 1
      ANG = FNF*SGN
      CSPNR = COS(ANG)
      CSPNI = SIN(ANG)
      IF (MOD(IFN,2).EQ.0) GO TO 190
      CSPNR = -CSPNR
      CSPNI = -CSPNI
  190 CONTINUE
C-----------------------------------------------------------------------
C     CS=COEFF OF THE J FUNCTION TO GET THE I FUNCTION. I(FNU,Z) IS
C     COMPUTED FROM EXP(I*FNU*HPI)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST
C     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY
C     CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS
C-----------------------------------------------------------------------
      CSR = SAR*CSGNI
      CSI = CAR*CSGNI
      IN = MOD(IFN,4) + 1
      C2R = CIPR(IN)
      C2I = CIPI(IN)
      STR = CSR*C2R + CSI*C2I
      CSI = -CSR*C2I + CSI*C2R
      CSR = STR
      ASC = BRY(1)
      IUF = 0
      KK = N
      KDFLG = 1
      IB = IB - 1
      IC = IB - 1
      DO 290 K=1,N
        FN = FNU + (KK-1)
C-----------------------------------------------------------------------
C     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
C     FUNCTION ABOVE
C-----------------------------------------------------------------------
        IF (N.GT.2) GO TO 175
  172   CONTINUE
        PHIDR = PHIR(J)
        PHIDI = PHII(J)
        ARGDR = ARGR(J)
        ARGDI = ARGI(J)
        ZET1DR = ZETA1R(J)
        ZET1DI = ZETA1I(J)
        ZET2DR = ZETA2R(J)
        ZET2DI = ZETA2I(J)
        ASUMDR = ASUMR(J)
        ASUMDI = ASUMI(J)
        BSUMDR = BSUMR(J)
        BSUMDI = BSUMI(J)
        J = 3 - J
        GO TO 210
  175   CONTINUE
        IF ((KK.EQ.N).AND.(IB.LT.N)) GO TO 210
        IF ((KK.EQ.IB).OR.(KK.EQ.IC)) GO TO 172
        CALL ZUNHJ(ZNR, ZNI, FN, 0, TOL, PHIDR, PHIDI, ARGDR,
     *   ARGDI, ZET1DR, ZET1DI, ZET2DR, ZET2DI, ASUMDR,
     *   ASUMDI, BSUMDR, BSUMDI)
  210   CONTINUE
        IF (KODE.EQ.1) GO TO 220
        STR = ZBR + ZET2DR
        STI = ZBI + ZET2DI
        RAST = FN/ZABS(STR,STI)
        STR = STR*RAST*RAST
        STI = -STI*RAST*RAST
        S1R = -ZET1DR + STR
        S1I = -ZET1DI + STI
        GO TO 230
  220   CONTINUE
        S1R = -ZET1DR + ZET2DR
        S1I = -ZET1DI + ZET2DI
  230   CONTINUE
C-----------------------------------------------------------------------
C     TEST FOR UNDERFLOW AND OVERFLOW
C-----------------------------------------------------------------------
        RS1 = S1R
        IF (ABS(RS1).GT.ELIM) GO TO 280
        IF (KDFLG.EQ.1) IFLAG = 2
        IF (ABS(RS1).LT.ALIM) GO TO 240
C-----------------------------------------------------------------------
C     REFINE  TEST AND SCALE
C-----------------------------------------------------------------------
        APHI = ZABS(PHIDR,PHIDI)
        AARG = ZABS(ARGDR,ARGDI)
        RS1 = RS1 + LOG(APHI) - 0.25D0*LOG(AARG) - AIC
        IF (ABS(RS1).GT.ELIM) GO TO 280
        IF (KDFLG.EQ.1) IFLAG = 1
        IF (RS1.LT.0.0D0) GO TO 240
        IF (KDFLG.EQ.1) IFLAG = 3
  240   CONTINUE
        CALL ZAIRY(ARGDR, ARGDI, 0, 2, AIR, AII, NAI, IDUM)
        CALL ZAIRY(ARGDR, ARGDI, 1, 2, DAIR, DAII, NDAI, IDUM)
        STR = DAIR*BSUMDR - DAII*BSUMDI
        STI = DAIR*BSUMDI + DAII*BSUMDR
        STR = STR + (AIR*ASUMDR-AII*ASUMDI)
        STI = STI + (AIR*ASUMDI+AII*ASUMDR)
        PTR = STR*PHIDR - STI*PHIDI
        PTI = STR*PHIDI + STI*PHIDR
        S2R = PTR*CSR - PTI*CSI
        S2I = PTR*CSI + PTI*CSR
        STR = EXP(S1R)*CSSR(IFLAG)
        S1R = STR*COS(S1I)
        S1I = STR*SIN(S1I)
        STR = S2R*S1R - S2I*S1I
        S2I = S2R*S1I + S2I*S1R
        S2R = STR
        IF (IFLAG.NE.1) GO TO 250
        CALL ZUCHK(S2R, S2I, NW, BRY(1), TOL)
        IF (NW.EQ.0) GO TO 250
        S2R = ZEROR
        S2I = ZEROI
  250   CONTINUE
        IF (YY.LE.0.0D0) S2I = -S2I
        CYR(KDFLG) = S2R
        CYI(KDFLG) = S2I
        C2R = S2R
        C2I = S2I
        S2R = S2R*CSRR(IFLAG)
        S2I = S2I*CSRR(IFLAG)
C-----------------------------------------------------------------------
C     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
C-----------------------------------------------------------------------
        S1R = YR(KK)
        S1I = YI(KK)
        IF (KODE.EQ.1) GO TO 270
        CALL ZS1S2(ZRR, ZRI, S1R, S1I, S2R, S2I, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  270   CONTINUE
        YR(KK) = S1R*CSPNR - S1I*CSPNI + S2R
        YI(KK) = S1R*CSPNI + S1I*CSPNR + S2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        STR = CSI
        CSI = -CSR
        CSR = STR
        IF (C2R.NE.0.0D0 .OR. C2I.NE.0.0D0) GO TO 255
        KDFLG = 1
        GO TO 290
  255   CONTINUE
        IF (KDFLG.EQ.2) GO TO 295
        KDFLG = 2
        GO TO 290
  280   CONTINUE
        IF (RS1.GT.0.0D0) GO TO 320
        S2R = ZEROR
        S2I = ZEROI
        GO TO 250
  290 CONTINUE
      K = N
  295 CONTINUE
      IL = N - K
      IF (IL.EQ.0) RETURN
C-----------------------------------------------------------------------
C     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
C     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
C     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
C-----------------------------------------------------------------------
      S1R = CYR(1)
      S1I = CYI(1)
      S2R = CYR(2)
      S2I = CYI(2)
      CSR = CSRR(IFLAG)
      ASCLE = BRY(IFLAG)
      FN = INU+IL
      DO 310 I=1,IL
        C2R = S2R
        C2I = S2I
        S2R = S1R + (FN+FNF)*(RZR*C2R-RZI*C2I)
        S2I = S1I + (FN+FNF)*(RZR*C2I+RZI*C2R)
        S1R = C2R
        S1I = C2I
        FN = FN - 1.0D0
        C2R = S2R*CSR
        C2I = S2I*CSR
        CKR = C2R
        CKI = C2I
        C1R = YR(KK)
        C1I = YI(KK)
        IF (KODE.EQ.1) GO TO 300
        CALL ZS1S2(ZRR, ZRI, C1R, C1I, C2R, C2I, NW, ASC, ALIM, IUF)
        NZ = NZ + NW
  300   CONTINUE
        YR(KK) = C1R*CSPNR - C1I*CSPNI + C2R
        YI(KK) = C1R*CSPNI + C1I*CSPNR + C2I
        KK = KK - 1
        CSPNR = -CSPNR
        CSPNI = -CSPNI
        IF (IFLAG.GE.3) GO TO 310
        C2R = ABS(CKR)
        C2I = ABS(CKI)
        C2M = MAX(C2R,C2I)
        IF (C2M.LE.ASCLE) GO TO 310
        IFLAG = IFLAG + 1
        ASCLE = BRY(IFLAG)
        S1R = S1R*CSR
        S1I = S1I*CSR
        S2R = CKR
        S2I = CKI
        S1R = S1R*CSSR(IFLAG)
        S1I = S1I*CSSR(IFLAG)
        S2R = S2R*CSSR(IFLAG)
        S2I = S2I*CSSR(IFLAG)
        CSR = CSRR(IFLAG)
  310 CONTINUE
      RETURN
  320 CONTINUE
      NZ = -1
      RETURN
      END
*DECK ZUOIK
      SUBROUTINE ZUOIK (ZR, ZI, FNU, KODE, IKFLG, N, YR, YI, NUF, TOL,
     +   ELIM, ALIM)
C***BEGIN PROLOGUE  ZUOIK
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESH, ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CUOIK-A, ZUOIK-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZUOIK COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC
C     EXPANSIONS FOR THE I AND K FUNCTIONS AND COMPARES THEM
C     (IN LOGARITHMIC FORM) TO ALIM AND ELIM FOR OVER AND UNDERFLOW
C     WHERE ALIM.LT.ELIM. IF THE MAGNITUDE, BASED ON THE LEADING
C     EXPONENTIAL, IS LESS THAN ALIM OR GREATER THAN -ALIM, THEN
C     THE RESULT IS ON SCALE. IF NOT, THEN A REFINED TEST USING OTHER
C     MULTIPLIERS (IN LOGARITHMIC FORM) IS MADE BASED ON ELIM. HERE
C     EXP(-ELIM)=SMALLEST MACHINE NUMBER*1.0E+3 AND EXP(-ALIM)=
C     EXP(-ELIM)/TOL
C
C     IKFLG=1 MEANS THE I SEQUENCE IS TESTED
C          =2 MEANS THE K SEQUENCE IS TESTED
C     NUF = 0 MEANS THE LAST MEMBER OF THE SEQUENCE IS ON SCALE
C         =-1 MEANS AN OVERFLOW WOULD OCCUR
C     IKFLG=1 AND NUF.GT.0 MEANS THE LAST NUF Y VALUES WERE SET TO ZERO
C             THE FIRST N-NUF VALUES MUST BE SET BY ANOTHER ROUTINE
C     IKFLG=2 AND NUF.EQ.N MEANS ALL Y VALUES WERE SET TO ZERO
C     IKFLG=2 AND 0.LT.NUF.LT.N NOT CONSIDERED. Y MUST BE SET BY
C             ANOTHER ROUTINE
C
C***SEE ALSO  ZBESH, ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZLOG, ZUCHK, ZUNHJ, ZUNIK
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C   930122  Added ZLOG to EXTERNAL statement.  (RWC)
C***END PROLOGUE  ZUOIK
C     COMPLEX ARG,ASUM,BSUM,CWRK,CZ,CZERO,PHI,SUM,Y,Z,ZB,ZETA1,ZETA2,ZN,
C    *ZR
      DOUBLE PRECISION AARG, AIC, ALIM, APHI, ARGI, ARGR, ASUMI, ASUMR,
     * ASCLE, AX, AY, BSUMI, BSUMR, CWRKI, CWRKR, CZI, CZR, ELIM, FNN,
     * FNU, GNN, GNU, PHII, PHIR, RCZ, STR, STI, SUMI, SUMR, TOL, YI,
     * YR, ZBI, ZBR, ZEROI, ZEROR, ZETA1I, ZETA1R, ZETA2I, ZETA2R, ZI,
     * ZNI, ZNR, ZR, ZRI, ZRR, D1MACH, ZABS
      INTEGER I, IDUM, IFORM, IKFLG, INIT, KODE, N, NN, NUF, NW
      DIMENSION YR(N), YI(N), CWRKR(16), CWRKI(16)
      EXTERNAL ZABS, ZLOG
      DATA ZEROR,ZEROI / 0.0D0, 0.0D0 /
      DATA AIC / 1.265512123484645396D+00 /
C***FIRST EXECUTABLE STATEMENT  ZUOIK
      NUF = 0
      NN = N
      ZRR = ZR
      ZRI = ZI
      IF (ZR.GE.0.0D0) GO TO 10
      ZRR = -ZR
      ZRI = -ZI
   10 CONTINUE
      ZBR = ZRR
      ZBI = ZRI
      AX = ABS(ZR)*1.7321D0
      AY = ABS(ZI)
      IFORM = 1
      IF (AY.GT.AX) IFORM = 2
      GNU = MAX(FNU,1.0D0)
      IF (IKFLG.EQ.1) GO TO 20
      FNN = NN
      GNN = FNU + FNN - 1.0D0
      GNU = MAX(GNN,FNN)
   20 CONTINUE
C-----------------------------------------------------------------------
C     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE
C     REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET
C     THE SIGN OF THE IMAGINARY PART CORRECT.
C-----------------------------------------------------------------------
      IF (IFORM.EQ.2) GO TO 30
      INIT = 0
      CALL ZUNIK(ZRR, ZRI, GNU, IKFLG, 1, TOL, INIT, PHIR, PHII,
     * ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
      CZR = -ZETA1R + ZETA2R
      CZI = -ZETA1I + ZETA2I
      GO TO 50
   30 CONTINUE
      ZNR = ZRI
      ZNI = -ZRR
      IF (ZI.GT.0.0D0) GO TO 40
      ZNR = -ZNR
   40 CONTINUE
      CALL ZUNHJ(ZNR, ZNI, GNU, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R,
     * ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
      CZR = -ZETA1R + ZETA2R
      CZI = -ZETA1I + ZETA2I
      AARG = ZABS(ARGR,ARGI)
   50 CONTINUE
      IF (KODE.EQ.1) GO TO 60
      CZR = CZR - ZBR
      CZI = CZI - ZBI
   60 CONTINUE
      IF (IKFLG.EQ.1) GO TO 70
      CZR = -CZR
      CZI = -CZI
   70 CONTINUE
      APHI = ZABS(PHIR,PHII)
      RCZ = CZR
C-----------------------------------------------------------------------
C     OVERFLOW TEST
C-----------------------------------------------------------------------
      IF (RCZ.GT.ELIM) GO TO 210
      IF (RCZ.LT.ALIM) GO TO 80
      RCZ = RCZ + LOG(APHI)
      IF (IFORM.EQ.2) RCZ = RCZ - 0.25D0*LOG(AARG) - AIC
      IF (RCZ.GT.ELIM) GO TO 210
      GO TO 130
   80 CONTINUE
C-----------------------------------------------------------------------
C     UNDERFLOW TEST
C-----------------------------------------------------------------------
      IF (RCZ.LT.(-ELIM)) GO TO 90
      IF (RCZ.GT.(-ALIM)) GO TO 130
      RCZ = RCZ + LOG(APHI)
      IF (IFORM.EQ.2) RCZ = RCZ - 0.25D0*LOG(AARG) - AIC
      IF (RCZ.GT.(-ELIM)) GO TO 110
   90 CONTINUE
      DO 100 I=1,NN
        YR(I) = ZEROR
        YI(I) = ZEROI
  100 CONTINUE
      NUF = NN
      RETURN
  110 CONTINUE
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      CALL ZLOG(PHIR, PHII, STR, STI, IDUM)
      CZR = CZR + STR
      CZI = CZI + STI
      IF (IFORM.EQ.1) GO TO 120
      CALL ZLOG(ARGR, ARGI, STR, STI, IDUM)
      CZR = CZR - 0.25D0*STR - AIC
      CZI = CZI - 0.25D0*STI
  120 CONTINUE
      AX = EXP(RCZ)/TOL
      AY = CZI
      CZR = AX*COS(AY)
      CZI = AX*SIN(AY)
      CALL ZUCHK(CZR, CZI, NW, ASCLE, TOL)
      IF (NW.NE.0) GO TO 90
  130 CONTINUE
      IF (IKFLG.EQ.2) RETURN
      IF (N.EQ.1) RETURN
C-----------------------------------------------------------------------
C     SET UNDERFLOWS ON I SEQUENCE
C-----------------------------------------------------------------------
  140 CONTINUE
      GNU = FNU + (NN-1)
      IF (IFORM.EQ.2) GO TO 150
      INIT = 0
      CALL ZUNIK(ZRR, ZRI, GNU, IKFLG, 1, TOL, INIT, PHIR, PHII,
     * ZETA1R, ZETA1I, ZETA2R, ZETA2I, SUMR, SUMI, CWRKR, CWRKI)
      CZR = -ZETA1R + ZETA2R
      CZI = -ZETA1I + ZETA2I
      GO TO 160
  150 CONTINUE
      CALL ZUNHJ(ZNR, ZNI, GNU, 1, TOL, PHIR, PHII, ARGR, ARGI, ZETA1R,
     * ZETA1I, ZETA2R, ZETA2I, ASUMR, ASUMI, BSUMR, BSUMI)
      CZR = -ZETA1R + ZETA2R
      CZI = -ZETA1I + ZETA2I
      AARG = ZABS(ARGR,ARGI)
  160 CONTINUE
      IF (KODE.EQ.1) GO TO 170
      CZR = CZR - ZBR
      CZI = CZI - ZBI
  170 CONTINUE
      APHI = ZABS(PHIR,PHII)
      RCZ = CZR
      IF (RCZ.LT.(-ELIM)) GO TO 180
      IF (RCZ.GT.(-ALIM)) RETURN
      RCZ = RCZ + LOG(APHI)
      IF (IFORM.EQ.2) RCZ = RCZ - 0.25D0*LOG(AARG) - AIC
      IF (RCZ.GT.(-ELIM)) GO TO 190
  180 CONTINUE
      YR(NN) = ZEROR
      YI(NN) = ZEROI
      NN = NN - 1
      NUF = NUF + 1
      IF (NN.EQ.0) RETURN
      GO TO 140
  190 CONTINUE
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      CALL ZLOG(PHIR, PHII, STR, STI, IDUM)
      CZR = CZR + STR
      CZI = CZI + STI
      IF (IFORM.EQ.1) GO TO 200
      CALL ZLOG(ARGR, ARGI, STR, STI, IDUM)
      CZR = CZR - 0.25D0*STR - AIC
      CZI = CZI - 0.25D0*STI
  200 CONTINUE
      AX = EXP(RCZ)/TOL
      AY = CZI
      CZR = AX*COS(AY)
      CZI = AX*SIN(AY)
      CALL ZUCHK(CZR, CZI, NW, ASCLE, TOL)
      IF (NW.NE.0) GO TO 180
      RETURN
  210 CONTINUE
      NUF = -1
      RETURN
      END
*DECK ZWRSK
      SUBROUTINE ZWRSK (ZRR, ZRI, FNU, KODE, N, YR, YI, NZ, CWR, CWI,
     +   TOL, ELIM, ALIM)
C***BEGIN PROLOGUE  ZWRSK
C***SUBSIDIARY
C***PURPOSE  Subsidiary to ZBESI and ZBESK
C***LIBRARY   SLATEC
C***TYPE      ALL (CWRSK-A, ZWRSK-A)
C***AUTHOR  Amos, D. E., (SNL)
C***DESCRIPTION
C
C     ZWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY
C     NORMALIZING THE I FUNCTION RATIOS FROM ZRATI BY THE WRONSKIAN
C
C***SEE ALSO  ZBESI, ZBESK
C***ROUTINES CALLED  D1MACH, ZABS, ZBKNU, ZRATI
C***REVISION HISTORY  (YYMMDD)
C   830501  DATE WRITTEN
C   910415  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZWRSK
C     COMPLEX CINU,CSCL,CT,CW,C1,C2,RCT,ST,Y,ZR
      DOUBLE PRECISION ACT, ACW, ALIM, ASCLE, CINUI, CINUR, CSCLR, CTI,
     * CTR, CWI, CWR, C1I, C1R, C2I, C2R, ELIM, FNU, PTI, PTR, RACT,
     * STI, STR, TOL, YI, YR, ZRI, ZRR, ZABS, D1MACH
      INTEGER I, KODE, N, NW, NZ
      DIMENSION YR(N), YI(N), CWR(2), CWI(2)
      EXTERNAL ZABS
C***FIRST EXECUTABLE STATEMENT  ZWRSK
C-----------------------------------------------------------------------
C     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
C     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
C     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
C-----------------------------------------------------------------------
C
      NZ = 0
      CALL ZBKNU(ZRR, ZRI, FNU, KODE, 2, CWR, CWI, NW, TOL, ELIM, ALIM)
      IF (NW.NE.0) GO TO 50
      CALL ZRATI(ZRR, ZRI, FNU, N, YR, YI, TOL)
C-----------------------------------------------------------------------
C     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
C     R(FNU+J-1,Z)=Y(J),  J=1,...,N
C-----------------------------------------------------------------------
      CINUR = 1.0D0
      CINUI = 0.0D0
      IF (KODE.EQ.1) GO TO 10
      CINUR = COS(ZRI)
      CINUI = SIN(ZRI)
   10 CONTINUE
C-----------------------------------------------------------------------
C     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
C     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
C     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
C     THE RESULT IS ON SCALE.
C-----------------------------------------------------------------------
      ACW = ZABS(CWR(2),CWI(2))
      ASCLE = 1.0D+3*D1MACH(1)/TOL
      CSCLR = 1.0D0
      IF (ACW.GT.ASCLE) GO TO 20
      CSCLR = 1.0D0/TOL
      GO TO 30
   20 CONTINUE
      ASCLE = 1.0D0/ASCLE
      IF (ACW.LT.ASCLE) GO TO 30
      CSCLR = TOL
   30 CONTINUE
      C1R = CWR(1)*CSCLR
      C1I = CWI(1)*CSCLR
      C2R = CWR(2)*CSCLR
      C2I = CWI(2)*CSCLR
      STR = YR(1)
      STI = YI(1)
C-----------------------------------------------------------------------
C     CINU=CINU*(CONJG(CT)/ABS(CT))*(1.0D0/ABS(CT) PREVENTS
C     UNDER- OR OVERFLOW PREMATURELY BY SQUARING ABS(CT)
C-----------------------------------------------------------------------
      PTR = STR*C1R - STI*C1I
      PTI = STR*C1I + STI*C1R
      PTR = PTR + C2R
      PTI = PTI + C2I
      CTR = ZRR*PTR - ZRI*PTI
      CTI = ZRR*PTI + ZRI*PTR
      ACT = ZABS(CTR,CTI)
      RACT = 1.0D0/ACT
      CTR = CTR*RACT
      CTI = -CTI*RACT
      PTR = CINUR*RACT
      PTI = CINUI*RACT
      CINUR = PTR*CTR - PTI*CTI
      CINUI = PTR*CTI + PTI*CTR
      YR(1) = CINUR*CSCLR
      YI(1) = CINUI*CSCLR
      IF (N.EQ.1) RETURN
      DO 40 I=2,N
        PTR = STR*CINUR - STI*CINUI
        CINUI = STR*CINUI + STI*CINUR
        CINUR = PTR
        STR = YR(I)
        STI = YI(I)
        YR(I) = CINUR*CSCLR
        YI(I) = CINUI*CSCLR
   40 CONTINUE
      RETURN
   50 CONTINUE
      NZ = -1
      IF(NW.EQ.(-2)) NZ=-2
      RETURN
      END
