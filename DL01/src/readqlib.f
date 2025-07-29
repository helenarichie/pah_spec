      SUBROUTINE READQLIB(ICASE,BOVERA,AEFF,WAVE,QABS,QEXT,QSCA)
      IMPLICIT NONE

!------------------------- readqlib_v19 ----------------------------
! arguments

      INTEGER ICASE
      REAL AEFF,BOVERA,WAVE
      REAL
     &   QABS(3),
     &   QEXT(3),
     &   QSCA(3)

! common/readqlibcom/ for testing purposes
      INTEGER RQTEST
      COMMON/READQLIBCOM/RQTEST

      CHARACTER*39 FLNAME
      COMMON/FLNAME_COM/FLNAME

! local variables
! NCTMAX = total number of compositions known to READQLIB
!          [currently 91]
! NCMAX  = maximum number of compositions that READQLIB can store
!          internally at one time (3 is a reasonable number)

      INTEGER NCMAX,NCTMAX,NRMAX,NWMAX
      PARAMETER(NWMAX=1128,NRMAX=168,NCMAX=3,NCTMAX=28+10*6+3)

      CHARACTER*1 C(0:9)
      CHARACTER*6 CHBA(-7:7)
      CHARACTER*17 CHPO(0:9)
      CHARACTER*14 CHPO2(0:9)
      CHARACTER*7 CHFE(0:5)
      CHARACTER*8 PATH
      CHARACTER*39 FILENAME(1:NCTMAX,-7:7)
      CHARACTER*37 FILENAME2(1:NCTMAX,-7:7)
      CHARACTER*20 TMPFILE
      INTEGER C1,CMAX,CRATE,IC,IDX,ISHAP,
     &   J1,J2,J3,J4,J5,J6,J7,J8,JBA,JC,JFE,JORI,JPO,JR,JS,JW,
     &   NCASE,NR,NW
      INTEGER 
     &   ICASE0(1:NCMAX),
     &   ICT(1:NCTMAX,-7:7),
     &   ISHAP0(1:NCMAX),
     &   NRAD(1:NCMAX),
     &   NWAV(1:NCMAX)
      REAL W00,W01,W10,W11,XR,XW
      REAL
     &   AXRAT(1:NCMAX),
     &   QTABS(1:3,0:NRMAX,0:NWMAX,1:NCMAX),
     &   QTEXT(1:3,0:NRMAX,0:NWMAX,1:NCMAX),
     &   QTSCA(1:3,0:NRMAX,0:NWMAX,1:NCMAX),
     &   RAD(0:NRMAX,1:NCMAX),
     &   WAV(0:NWMAX,1:NCMAX)

! 169 sizes (40 per decade, 3.162e-4um to 5.012um)
! 1129 wavelengths (200 per decade, .09120um to 3.981e4um)
! array size: 3*(NRMAX+1)*(NWMAX+1)*NCMAX*4 bytes
!            =3*   169   *  1129   *  3  *4 bytes = 6.87 Mbytes

      DATA C/'0','1','2','3','4','5','6','7','8','9'/

      SAVE FILENAME,PATH
      SAVE ICASE0,ICT,ISHAP0,NCASE,NRAD,NWAV
      SAVE QTABS,QTEXT,QTSCA,RAD,WAV

!----------------- readqlib_v19 ------------------------------------
! subroutine READQLIB
! given:
!       ICASE = 1 for Draine 2003 astrosilicate
!               2 for graphite-PAH^+,a_t=50A
!               3 for graphite-PAH^0,a_t=50A
!               4 for HD2014 silicate, no Fe incl
!               5 for HD2014 sil + Fe incl, fill=0.02
!               6 for HD2014 sil + Fe incl, fill=0.05
!               7 for HD2014 sil + Fe incl, fill=0.10
!               8 for HD2014 sil + Fe incl, fill=0.20
!               9 for DH2013 Fe only
!              10 for HD2014 graphite-PAH^+ [reduced sigma_dc(E||c)],a_t=20A
!              11 for HD2014 graphite-PAH^0 [reduced sigma_dc(E||c)],a_t=20A
!              12 for HD2014 graphite-PAH^+ [reduced sigma_dc(E||c)],a_t=50A
!              13 for HD2014 graphite-PAH^0 [reduced sigma_dc(E||c)],a_t=50A
!              14 for D16 graphite-PAH+, random c axis, a_t=20A
!              15 for D16 graphite-PAH0, random c axis, a_t=20A
!              16 for D16 graphite-PAH+, random c axis, a_t=50A
!              17 for D16 graphite-PAH0, random c axis, a_t=50A
!              18 for D16 graphite-PAH+, c || symm axis a, a_t=20A
!              19 for D16 graphite-PAH0, c || symm axis a, a_t=20A
!              20 for D16 graphite-PAH+, c || symm axis a, a_t=50A
!              21 for D16 turbost. graphite-PAH0, c || symm axis a, a_t=50A
!              22 for D16 turbost. graphite-PAH+, MG1 EMT (matrix=E||c), a_t=20A
!              23 for D16 turbost. graphite-PAH0, MG1 EMT (matrix=E||c), a_t=20A
!              24 for D16 turbost. graphite-PAH+, MG1 EMT (matrix=E||c), a_t=50A
!              25 for D16 turbost. graphite-PAH0, MG1 EMT (matrix=E||c), a_t=50A
!              26 for D16 graphite, random c axis
!              27 for D16 graphite, c || a
!              28 for D16 turbostratic graphite, EMT (MG1: matrix=E||c)
!           29-88 for DH21 Astrodust
!                 ICASE=29+IDX
!                 IDX= JPO*6 + JFE
!                 JPO=0-9 for porosity= 0, 0.1, 0.2, ..., 0.9  ]10 cases]
!                 JFE=0-5 for f_Fe=0, 0.10, 0.20, ... 0.50   [ 6 cases]
!              89 for D16 graphite, MG2 EMT (MG2: MG with matrix=Eperpc)
!              90 for D16 graphite, MG2 EMT-PAH+, a_t=50A
!              91 for D16 graphite, MG2 EMT-PAH0, a_t=50A
!      examples:

!              29 for porosity = 0  , f_Fe = 0

!              30 for porosity = 0  , f_Fe = 0.1

!              35 for porosity = 0.1, f_Fe = 0
!              41                0.2         0
!              47                0.3         0
!              53                0.4         0
!              59                0.5         0
!              65                0.6         0
!              71                0.7         0
!              77                0.8         0
!              83                0.9         0
      
!              88 for porosity = 0.9, f_Fe = 0.5

!       BOVERA = b/a
!               allowed values:
!                  0.333, 0.400, 0.500, 0.556, 0.625, 0.714, 0.833, 1.000,
!                  1.200, 1.400, 1.600, 1.800, 2.000, 2.500, 3.000

!       AEFF  = aeff = (a*b*b)^{1/3} (um)
!       WAVE  = wavelength (um)
! returns
!       QABS(1-3) = Cabs/pi*aeff^2 for ori=1,2,3
!       QEXT(1-3) = Cext/pi*aeff^2 for ori=1,2,3
!       QSCA(1-3) = Csca/pi*aeff^2 for ori=1,2,3

! Let a = symmetry axis
! iori=1: k  ||  a  (and E perp a)
!      2: k perp a, E  ||  a  (TM mode)
!      3: k perp a, E perp a  (TE mode)

! Q values are obtained by table lookup/interpolation
!          with extrapolation if size or wavelength are
!          out of bounds
! stored tables generated by program /u/draine/work/opt/exe/callmakeqlib
!
! B.T. Draine, Princeton University, 2013.07.09
! history
! 13.07.09 (BTD) minor edits for first version.
! 13.07.12 (BTD) switch positions of QABS and QEXT in arg list
!                they are now in alphabetical order: QABS,QEXT,QSCA
!                modify how cases are stored
! 13.09.30 (BTD) changed directory from opt/dat to opt/qdat
! 14.02.12 (BTD) v2 
!                * interim change to accept ICASE,ISHAP as arguments.
!                * change to new file names
! 14.02.12 (BTD) v3
!                * change to correctly handle two-dimensional grid of
!                  models, parameterized by
!                  ICASE for composition
!                  ISHAP for shape
!                * multiple changes made to accomodate new character
! 14.04.24 (BTD) * added comments describing iori=1,2,3
! 14.07.22 (BTD) v4
!                * added support for case10 and case11
! 14.07.24 (BTD) * added support for case12 and case13
! 14.09.03 (BTD) v5
!                * added support for ICASE=14, 15, 16, 17, 18
! 16.02.08 (BTD) v6
!                * added support for ICASE=19, 20, 21, 22
! 16.02.13 (BTD) v7
!                * added support for ICASE=23, 24, 25, 26
!                                          27, 28, 29, 30
! 16.02.19 (BTD) * added STATUS='OLD' to OPEN statement
! 16.09.29 (BTD) * added trap to warn if used outside
!                  range of tabulated sizes or wavelengths
! 17.04.06 (BTD) v8 modified to extend ICASE to 38
! 17.05.19 (BTD) v9 modified to extend shape to
!                   JC=6 : b/a=2.5
!                      7 : b/a=3
!                   and add loop to initialize file names to
!                   'not yet defined' until specifically redefined
! 17.07.19 (BTD)   modified to extend ICASE to 40
! 17.09.22 (BTD) v10
!                modified to use new set of cases
!                extend to have 15 standard axial ratios from 0.333 to 3
! 17.12.21 (BTD) modified to treat cases ICASE=41 and 42
! 18.01.25 (BTD) v12
!                * extend to treat ICASE up to 80
!                * renumber ICASE=41,42 -> 79,80
!                * rename DH2017 -> DH2018
!                * replaced executable definition of FILENAME
!                  by DATA statements
! 18.01.26 (BTD) v13
!                * modified to read gzipped qlib files
!                  using system call to gunzip the desired file
!                  into a tempfile that is deleted when finished
! 18.06.26 (BTD) v14
!                * modified to return filename through 
!                  COMMON/FLNAME_COM/FLNAME to support reading of
!                  qlib_index file
! 18.06.27 (BTD) v15
!                * major rewrite to switch to new parameterization
!                  based on porosity rather than alphasil
!                  and on fraction of solid volume occupied by Fe
!                  rather than fraction of Fe in metallic form
! 18.08.26 (BTD) * changed FLNAME from CHARACTER*38 to CHARACTER*39
! 19.01.30 (BTD) v17
!                * modified filename coverage for ICASE=29-568
!                  JF8 now runs from 4 to 12
!                  f8 runs from 0.4 to 1.2
!                * support for 
!                  ICASE=569: D16 graphite, MG2 EMT
!                  ICASE=570: D16 graphite, MG2 EMT, PAH+, a_t=50A
!                  ICASE=571: D16 graphite, MG2 EMT, PAH0, a_t=50A
! 19.07.25 (BTD) v18
!                * modified filename coverage for ICASE=29-88,89,90,91
! 19.09.06 (BTD) * add ICASE0 and ISHAP0 to SAVE statement
! 20.06.09 (BTD) * changed filenames
!                  qlib_gra_d16emt2pahib > qlib_gra_d16emt2_pahib
!                  qlib_gra_d16emt2pahnb > qlib_gra_d16emt2_pahnb
! 20.12.03 (BTD) v19
!                * change file names to DH21Ad...
!                * for DH21Ad cases: 
!                  on first call (when reading in original qlib file)
!                  write out files 
!                  DH21_aeff with list of radii
!                  DH21_wave with list of wavelengths
!                  q_DH21Ad_Pxxx_Fexxx_bovera.gz
!                           where Pxxx =  P0.00, P0.10, P0.20...
!                                 Fexxx= Fe0.00, Fe0.10, ...  
!                                 bovera=0.500, 0.625, 1.400, 2.000, ...
!                           and file is gzipped
! 21.02.03 (BTD) * corrected calculation of NCTMAX
!                  (was larger than needed)
! end history
!=======================================================================
      DATA NCASE/0/
      DATA (FILENAME(1,JBA),JBA=-7,7)/
     &  'qlib_silicate_0.333.gz',
     &'  qlib_silicate_0.400.gz',
     &  'qlib_silicate_0.500.gz',
     &  'qlib_silicate_0.556.gz',
     &  'qlib_silicate_0.625.gz',
     &  'qlib_silicate_0.714.gz',
     &  'qlib_silicate_0.833.gz',
     &  'qlib_silicate_1.000.gz',
     &  'qlib_silicate_1.200.gz',
     &  'qlib_silicate_1.400.gz',
     &  'qlib_silicate_1.600.gz',
     &  'qlib_silicate_1.800.gz',
     &  'qlib_silicate_2.000.gz',
     &  'qlib_silicate_2.500.gz',
     &  'qlib_silicate_3.000.gz'/

      DATA (FILENAME(2,JBA),JBA=-7,7)/
     &   'qlib_gra-pahi_0.333.gz',
     &   'qlib_gra-pahi_0.400.gz',
     &   'qlib_gra-pahi_0.500.gz',
     &   'qlib_gra-pahi_0.556.gz',
     &   'qlib_gra-pahi_0.625.gz',
     &   'qlib_gra-pahi_0.714.gz',
     &   'qlib_gra-pahi_0.833.gz',
     &   'qlib_gra-pahi_1.000.gz',
     &   'qlib_gra-pahi_1.200.gz',
     &   'qlib_gra-pahi_1.400.gz',
     &   'qlib_gra-pahi_1.600.gz',
     &   'qlib_gra-pahi_1.800.gz',
     &   'qlib_gra-pahi_2.000.gz',
     &   'qlib_gra-pahi_2.500.gz',
     &   'qlib_gra-pahi_3.000.gz'/

      DATA (FILENAME(3,JBA),JBA=-7,7)/
     &   'qlib_gra-pahn_0.333.gz',
     &   'qlib_gra-pahn_0.400.gz',
     &   'qlib_gra-pahn_0.500.gz',
     &   'qlib_gra-pahn_0.556.gz',
     &   'qlib_gra-pahn_0.625.gz',
     &   'qlib_gra-pahn_0.714.gz',
     &   'qlib_gra-pahn_0.833.gz',
     &   'qlib_gra-pahn_1.000.gz',
     &   'qlib_gra-pahn_1.200.gz',
     &   'qlib_gra-pahn_1.400.gz',
     &   'qlib_gra-pahn_1.600.gz',
     &   'qlib_gra-pahn_1.800.gz',
     &   'qlib_gra-pahn_2.000.gz',
     &   'qlib_gra-pahn_2.500.gz',
     &   'qlib_gra-pahn_3.000.gz'/

      DATA (FILENAME(4,JBA),JBA=-7,7)/
     &   'qlib_sil_HD14_0.333.gz',
     &   'qlib_sil_HD14_0.400.gz',
     &   'qlib_sil_HD14_0.500.gz',
     &   'qlib_sil_HD14_0.556.gz',
     &   'qlib_sil_HD14_0.625.gz',
     &   'qlib_sil_HD14_0.714.gz',
     &   'qlib_sil_HD14_0.833.gz',
     &   'qlib_sil_HD14_1.000.gz',
     &   'qlib_sil_HD14_1.200.gz',
     &   'qlib_sil_HD14_1.400.gz',
     &   'qlib_sil_HD14_1.600.gz',
     &   'qlib_sil_HD14_1.800.gz',
     &   'qlib_sil_HD14_2.000.gz',
     &   'qlib_sil_HD14_2.500.gz',
     &   'qlib_sil_HD14_3.000.gz'/

      DATA (FILENAME(5,JBA),JBA=-7,7)/
     &   'qlib_sil_HD14_Fe0.02_0.333.gz',
     &   'qlib_sil_HD14_Fe0.02_0.400.gz',
     &   'qlib_sil_HD14_Fe0.02_0.500.gz',
     &   'qlib_sil_HD14_Fe0.02_0.556.gz',
     &   'qlib_sil_HD14_Fe0.02_0.625.gz',
     &   'qlib_sil_HD14_Fe0.02_0.714.gz',
     &   'qlib_sil_HD14_Fe0.02_0.833.gz',
     &   'qlib_sil_HD14_Fe0.02_1.000.gz',
     &   'qlib_sil_HD14_Fe0.02_1.200.gz',
     &   'qlib_sil_HD14_Fe0.02_1.400.gz',
     &   'qlib_sil_HD14_Fe0.02_1.600.gz',
     &   'qlib_sil_HD14_Fe0.02_1.800.gz',
     &   'qlib_sil_HD14_Fe0.02_2.000.gz',
     &   'qlib_sil_HD14_Fe0.02_2.500.gz',
     &   'qlib_sil_HD14_Fe0.02_3.000.gz'/
 
       DATA (FILENAME(6,JBA),JBA=-7,7)/
     &   'qlib_sil_HD14_Fe0.05_0.333.gz',
     &   'qlib_sil_HD14_Fe0.05_0.400.gz',
     &   'qlib_sil_HD14_Fe0.05_0.500.gz',
     &   'qlib_sil_HD14_Fe0.05_0.556.gz',
     &   'qlib_sil_HD14_Fe0.05_0.625.gz',
     &   'qlib_sil_HD14_Fe0.05_0.714.gz',
     &   'qlib_sil_HD14_Fe0.05_0.833.gz',
     &   'qlib_sil_HD14_Fe0.05_1.000.gz',
     &   'qlib_sil_HD14_Fe0.05_1.200.gz',
     &   'qlib_sil_HD14_Fe0.05_1.400.gz',
     &   'qlib_sil_HD14_Fe0.05_1.600.gz',
     &   'qlib_sil_HD14_Fe0.05_1.800.gz',
     &   'qlib_sil_HD14_Fe0.05_2.000.gz',
     &   'qlib_sil_HD14_Fe0.05_2.500.gz',
     &   'qlib_sil_HD14_Fe0.05_3.000.gz'/

       DATA (FILENAME(7,JBA),JBA=-7,7)/
     &   'qlib_sil_HD14_Fe0.10_0.333.gz',
     &   'qlib_sil_HD14_Fe0.10_0.400.gz',
     &   'qlib_sil_HD14_Fe0.10_0.500.gz',
     &   'qlib_sil_HD14_Fe0.10_0.556.gz',
     &   'qlib_sil_HD14_Fe0.10_0.625.gz',
     &   'qlib_sil_HD14_Fe0.10_0.714.gz',
     &   'qlib_sil_HD14_Fe0.10_0.833.gz',
     &   'qlib_sil_HD14_Fe0.10_1.000.gz',
     &   'qlib_sil_HD14_Fe0.10_1.200.gz',
     &   'qlib_sil_HD14_Fe0.10_1.400.gz',
     &   'qlib_sil_HD14_Fe0.10_1.600.gz',
     &   'qlib_sil_HD14_Fe0.10_1.800.gz',
     &   'qlib_sil_HD14_Fe0.10_2.000.gz',
     &   'qlib_sil_HD14_Fe0.10_2.500.gz',
     &   'qlib_sil_HD14_Fe0.10_3.000.gz'/

       DATA (FILENAME(8,JBA),JBA=-7,7)/
     &   'qlib_sil_HD14_Fe0.20_0.333.gz',
     &   'qlib_sil_HD14_Fe0.20_0.400.gz',
     &   'qlib_sil_HD14_Fe0.20_0.500.gz',
     &   'qlib_sil_HD14_Fe0.20_0.556.gz',
     &   'qlib_sil_HD14_Fe0.20_0.625.gz',
     &   'qlib_sil_HD14_Fe0.20_0.714.gz',
     &   'qlib_sil_HD14_Fe0.20_0.833.gz',
     &   'qlib_sil_HD14_Fe0.20_1.000.gz',
     &   'qlib_sil_HD14_Fe0.20_1.200.gz',
     &   'qlib_sil_HD14_Fe0.20_1.400.gz',
     &   'qlib_sil_HD14_Fe0.20_1.600.gz',
     &   'qlib_sil_HD14_Fe0.20_1.800.gz',
     &   'qlib_sil_HD14_Fe0.20_2.000.gz',
     &   'qlib_sil_HD14_Fe0.20_2.500.gz',
     &   'qlib_sil_HD14_Fe0.20_3.000.gz'/

      DATA (FILENAME(9,JBA),JBA=-7,7)/
     &   'qlib_DH13_Fe_0.333.gz',
     &   'qlib_DH13_Fe_0.400.gz',
     &   'qlib_DH13_Fe_0.500.gz',
     &   'qlib_DH13_Fe_0.556.gz',
     &   'qlib_DH13_Fe_0.625.gz',
     &   'qlib_DH13_Fe_0.714.gz',
     &   'qlib_DH13_Fe_0.833.gz',
     &   'qlib_DH13_Fe_1.000.gz',
     &   'qlib_DH13_Fe_1.200.gz',
     &   'qlib_DH13_Fe_1.400.gz',
     &   'qlib_DH13_Fe_1.600.gz',
     &   'qlib_DH13_Fe_1.800.gz',
     &   'qlib_DH13_Fe_2.000.gz',
     &   'qlib_DH13_Fe_2.500.gz',
     &   'qlib_DH13_Fe_3.000.gz'/

      DATA (FILENAME(10,JBA),JBA=-7,7)/
     &   'qlib_gra-pahi_HD14_0.333.gz',
     &   'qlib_gra-pahi_HD14_0.400.gz',
     &   'qlib_gra-pahi_HD14_0.500.gz',
     &   'qlib_gra-pahi_HD14_0.556.gz',
     &   'qlib_gra-pahi_HD14_0.625.gz',
     &   'qlib_gra-pahi_HD14_0.714.gz',
     &   'qlib_gra-pahi_HD14_0.833.gz',
     &   'qlib_gra-pahi_HD14_1.000.gz',
     &   'qlib_gra-pahi_HD14_1.200.gz',
     &   'qlib_gra-pahi_HD14_1.400.gz',
     &   'qlib_gra-pahi_HD14_1.600.gz',
     &   'qlib_gra-pahi_HD14_1.800.gz',
     &   'qlib_gra-pahi_HD14_2.000.gz',
     &   'qlib_gra-pahi_HD14_2.500.gz',
     &   'qlib_gra-pahi_HD14_3.000.gz'/

      DATA (FILENAME(11,JBA),JBA=-7,7)/
     &   'qlib_gra-pahn_HD14_0.333.gz',
     &   'qlib_gra-pahn_HD14_0.400.gz',
     &   'qlib_gra-pahn_HD14_0.500.gz',
     &   'qlib_gra-pahn_HD14_0.556.gz',
     &   'qlib_gra-pahn_HD14_0.625.gz',
     &   'qlib_gra-pahn_HD14_0.714.gz',
     &   'qlib_gra-pahn_HD14_0.833.gz',
     &   'qlib_gra-pahn_HD14_1.000.gz',
     &   'qlib_gra-pahn_HD14_1.200.gz',
     &   'qlib_gra-pahn_HD14_1.400.gz',
     &   'qlib_gra-pahn_HD14_1.600.gz',
     &   'qlib_gra-pahn_HD14_1.800.gz',
     &   'qlib_gra-pahn_HD14_2.000.gz',
     &   'qlib_gra-pahn_HD14_2.500.gz',
     &   'qlib_gra-pahn_HD14_3.000.gz'/

 
      DATA (FILENAME(12,JBA),JBA=-7,7)/
     &   'qlib_gra-pahi_HD14b_0.333.gz',
     &   'qlib_gra-pahi_HD14b_0.400.gz',
     &   'qlib_gra-pahi_HD14b_0.500.gz',
     &   'qlib_gra-pahi_HD14b_0.556.gz',
     &   'qlib_gra-pahi_HD14b_0.625.gz',
     &   'qlib_gra-pahi_HD14b_0.714.gz',
     &   'qlib_gra-pahi_HD14b_0.833.gz',
     &   'qlib_gra-pahi_HD14b_1.000.gz',
     &   'qlib_gra-pahi_HD14b_1.200.gz',
     &   'qlib_gra-pahi_HD14b_1.400.gz',
     &   'qlib_gra-pahi_HD14b_1.600.gz',
     &   'qlib_gra-pahi_HD14b_1.800.gz',
     &   'qlib_gra-pahi_HD14b_2.000.gz',
     &   'qlib_gra-pahi_HD14b_2.500.gz',
     &   'qlib_gra-pahi_HD14b_3.000.gz'/

      DATA (FILENAME(13,JBA),JBA=-7,7)/
     &   'qlib_gra-pahn_HD14b_0.333.gz',
     &   'qlib_gra-pahn_HD14b_0.400.gz',
     &   'qlib_gra-pahn_HD14b_0.500.gz',
     &   'qlib_gra-pahn_HD14b_0.556.gz',
     &   'qlib_gra-pahn_HD14b_0.625.gz',
     &   'qlib_gra-pahn_HD14b_0.714.gz',
     &   'qlib_gra-pahn_HD14b_0.833.gz',
     &   'qlib_gra-pahn_HD14b_1.000.gz',
     &   'qlib_gra-pahn_HD14b_1.200.gz',
     &   'qlib_gra-pahn_HD14b_1.400.gz',
     &   'qlib_gra-pahn_HD14b_1.600.gz',
     &   'qlib_gra-pahn_HD14b_1.800.gz',
     &   'qlib_gra-pahn_HD14b_2.000.gz',
     &   'qlib_gra-pahn_HD14b_2.500.gz',
     &   'qlib_gra-pahn_HD14b_3.000.gz'/

! D16 graphite, PAH+, a_t=20A, f_gmin=0.01
!               randomly-oriented c axis
      DATA (FILENAME(14,JBA),JBA=-7,7)/
     &   'qlib_gra_D16rc_pahi_0.333.gz',
     &   'qlib_gra_D16rc_pahi_0.400.gz',
     &   'qlib_gra_D16rc_pahi_0.500.gz',
     &   'qlib_gra_D16rc_pahi_0.556.gz',
     &   'qlib_gra_D16rc_pahi_0.625.gz',
     &   'qlib_gra_D16rc_pahi_0.714.gz',
     &   'qlib_gra_D16rc_pahi_0.833.gz',
     &   'qlib_gra_D16rc_pahi_1.000.gz',
     &   'qlib_gra_D16rc_pahi_1.200.gz',
     &   'qlib_gra_D16rc_pahi_1.400.gz',
     &   'qlib_gra_D16rc_pahi_1.600.gz',
     &   'qlib_gra_D16rc_pahi_1.800.gz',
     &   'qlib_gra_D16rc_pahi_2.000.gz',
     &   'qlib_gra_D16rc_pahi_2.500.gz',
     &   'qlib_gra_D16rc_pahi_3.000.gz'/

! D16 graphite, PAH0, a_t=20A, f_gmin=0.01
!               randomly-oriented c axis
      DATA (FILENAME(15,JBA),JBA=-7,7)/
     &   'qlib_gra_D16rc_pahn_0.333.gz',
     &   'qlib_gra_D16rc_pahn_0.400.gz',
     &   'qlib_gra_D16rc_pahn_0.500.gz',
     &   'qlib_gra_D16rc_pahn_0.556.gz',
     &   'qlib_gra_D16rc_pahn_0.625.gz',
     &   'qlib_gra_D16rc_pahn_0.714.gz',
     &   'qlib_gra_D16rc_pahn_0.833.gz',
     &   'qlib_gra_D16rc_pahn_1.000.gz',
     &   'qlib_gra_D16rc_pahn_1.200.gz',
     &   'qlib_gra_D16rc_pahn_1.400.gz',
     &   'qlib_gra_D16rc_pahn_1.600.gz',
     &   'qlib_gra_D16rc_pahn_1.800.gz',
     &   'qlib_gra_D16rc_pahn_2.000.gz',
     &   'qlib_gra_D16rc_pahn_2.500.gz',
     &   'qlib_gra_D16rc_pahn_3.000.gz'/

! D16 graphite, PAH+, a_t=50A, f_gmin=0.01
!               randomly-oriented c axis
      DATA (FILENAME(16,JBA),JBA=-7,7)/
     &   'qlib_gra_D16rc_pahib_0.333.gz',
     &   'qlib_gra_D16rc_pahib_0.400.gz',
     &   'qlib_gra_D16rc_pahib_0.500.gz',
     &   'qlib_gra_D16rc_pahib_0.556.gz',
     &   'qlib_gra_D16rc_pahib_0.625.gz',
     &   'qlib_gra_D16rc_pahib_0.714.gz',
     &   'qlib_gra_D16rc_pahib_0.833.gz',
     &   'qlib_gra_D16rc_pahib_1.000.gz',
     &   'qlib_gra_D16rc_pahib_1.200.gz',
     &   'qlib_gra_D16rc_pahib_1.400.gz',
     &   'qlib_gra_D16rc_pahib_1.600.gz',
     &   'qlib_gra_D16rc_pahib_1.800.gz',
     &   'qlib_gra_D16rc_pahib_2.000.gz',
     &   'qlib_gra_D16rc_pahib_2.500.gz',
     &   'qlib_gra_D16rc_pahib_3.000.gz'/

! D16 graphite, PAH0, a_t=20A, f_gmin=0.01
!               randomly-oriented c axis
      DATA (FILENAME(17,JBA),JBA=-7,7)/
     &   'qlib_gra_D16rc_pahnb_0.333.gz',
     &   'qlib_gra_D16rc_pahnb_0.400.gz',
     &   'qlib_gra_D16rc_pahnb_0.500.gz',
     &   'qlib_gra_D16rc_pahnb_0.556.gz',
     &   'qlib_gra_D16rc_pahnb_0.625.gz',
     &   'qlib_gra_D16rc_pahnb_0.714.gz',
     &   'qlib_gra_D16rc_pahnb_0.833.gz',
     &   'qlib_gra_D16rc_pahnb_1.000.gz',
     &   'qlib_gra_D16rc_pahnb_1.200.gz',
     &   'qlib_gra_D16rc_pahnb_1.400.gz',
     &   'qlib_gra_D16rc_pahnb_1.600.gz',
     &   'qlib_gra_D16rc_pahnb_1.800.gz',
     &   'qlib_gra_D16rc_pahnb_2.000.gz',
     &   'qlib_gra_D16rc_pahnb_2.500.gz',
     &   'qlib_gra_D16rc_pahnb_3.000.gz'/

! D16 graphite, PAH+, a_t=20A, f_gmin=0.01
!               c axis || symm axis a
      DATA (FILENAME(18,JBA),JBA=-7,7)/
     &   'qlib_gra_D16ca_pahi_0.333.gz',
     &   'qlib_gra_D16ca_pahi_0.400.gz',
     &   'qlib_gra_D16ca_pahi_0.500.gz',
     &   'qlib_gra_D16ca_pahi_0.556.gz',
     &   'qlib_gra_D16ca_pahi_0.625.gz',
     &   'qlib_gra_D16ca_pahi_0.714.gz',
     &   'qlib_gra_D16ca_pahi_0.833.gz',
     &   'qlib_gra_D16ca_pahi_1.000.gz',
     &   'qlib_gra_D16ca_pahi_1.200.gz',
     &   'qlib_gra_D16ca_pahi_1.400.gz',
     &   'qlib_gra_D16ca_pahi_1.600.gz',
     &   'qlib_gra_D16ca_pahi_1.800.gz',
     &   'qlib_gra_D16ca_pahi_2.000.gz',
     &   'qlib_gra_D16ca_pahi_2.500.gz',
     &   'qlib_gra_D16ca_pahi_3.000.gz'/

! D16 graphite, PAH0, a_t=20A, f_gmin=0.01
!               c axis || symm axis a
      DATA (FILENAME(19,JBA),JBA=-7,7)/
     &   'qlib_gra_D16ca_pahn_0.333.gz',
     &   'qlib_gra_D16ca_pahn_0.400.gz',
     &   'qlib_gra_D16ca_pahn_0.500.gz',
     &   'qlib_gra_D16ca_pahn_0.556.gz',
     &   'qlib_gra_D16ca_pahn_0.625.gz',
     &   'qlib_gra_D16ca_pahn_0.714.gz',
     &   'qlib_gra_D16ca_pahn_0.833.gz',
     &   'qlib_gra_D16ca_pahn_1.000.gz',
     &   'qlib_gra_D16ca_pahn_1.200.gz',
     &   'qlib_gra_D16ca_pahn_1.400.gz',
     &   'qlib_gra_D16ca_pahn_1.600.gz',
     &   'qlib_gra_D16ca_pahn_1.800.gz',
     &   'qlib_gra_D16ca_pahn_2.000.gz',
     &   'qlib_gra_D16ca_pahn_2.500.gz',
     &   'qlib_gra_D16ca_pahn_3.000.gz'/

! D16 graphite, PAH+, a_t=50A, f_gmin=0.01
!               c axis || symm axis a
      DATA (FILENAME(20,JBA),JBA=-7,7)/
     &   'qlib_gra_D16ca_pahib_0.333.gz',
     &   'qlib_gra_D16ca_pahib_0.400.gz',
     &   'qlib_gra_D16ca_pahib_0.500.gz',
     &   'qlib_gra_D16ca_pahib_0.556.gz',
     &   'qlib_gra_D16ca_pahib_0.625.gz',
     &   'qlib_gra_D16ca_pahib_0.714.gz',
     &   'qlib_gra_D16ca_pahib_0.833.gz',
     &   'qlib_gra_D16ca_pahib_1.000.gz',
     &   'qlib_gra_D16ca_pahib_1.200.gz',
     &   'qlib_gra_D16ca_pahib_1.400.gz',
     &   'qlib_gra_D16ca_pahib_1.600.gz',
     &   'qlib_gra_D16ca_pahib_1.800.gz',
     &   'qlib_gra_D16ca_pahib_2.000.gz',
     &   'qlib_gra_D16ca_pahib_2.500.gz',
     &   'qlib_gra_D16ca_pahib_3.000.gz'/

! D16 graphite, PAH0, a_t=50A, f_gmin=0.01
!               c axis || symm axis a
      DATA (FILENAME(21,JBA),JBA=-7,7)/
     &   'qlib_gra_D16ca_pahnb_0.333.gz',
     &   'qlib_gra_D16ca_pahnb_0.400.gz',
     &   'qlib_gra_D16ca_pahnb_0.500.gz',
     &   'qlib_gra_D16ca_pahnb_0.556.gz',
     &   'qlib_gra_D16ca_pahnb_0.625.gz',
     &   'qlib_gra_D16ca_pahnb_0.714.gz',
     &   'qlib_gra_D16ca_pahnb_0.833.gz',
     &   'qlib_gra_D16ca_pahnb_1.000.gz',
     &   'qlib_gra_D16ca_pahnb_1.200.gz',
     &   'qlib_gra_D16ca_pahnb_1.400.gz',
     &   'qlib_gra_D16ca_pahnb_1.600.gz',
     &   'qlib_gra_D16ca_pahnb_1.800.gz',
     &   'qlib_gra_D16ca_pahnb_2.000.gz',
     &   'qlib_gra_D16ca_pahnb_2.500.gz',
     &   'qlib_gra_D16ca_pahnb_3.000.gz'/

! D16 graphite, PAH+, a_t=20A, f_gmin=0.01
!               EMT used for graphite
      DATA (FILENAME(22,JBA),JBA=-7,7)/
     &   'qlib_gra_D16emt_pahi_0.333.gz',
     &   'qlib_gra_D16emt_pahi_0.400.gz',
     &   'qlib_gra_D16emt_pahi_0.500.gz',
     &   'qlib_gra_D16emt_pahi_0.556.gz',
     &   'qlib_gra_D16emt_pahi_0.625.gz',
     &   'qlib_gra_D16emt_pahi_0.714.gz',
     &   'qlib_gra_D16emt_pahi_0.833.gz',
     &   'qlib_gra_D16emt_pahi_1.000.gz',
     &   'qlib_gra_D16emt_pahi_1.200.gz',
     &   'qlib_gra_D16emt_pahi_1.400.gz',
     &   'qlib_gra_D16emt_pahi_1.600.gz',
     &   'qlib_gra_D16emt_pahi_1.800.gz',
     &   'qlib_gra_D16emt_pahi_2.000.gz',
     &   'qlib_gra_D16emt_pahi_2.500.gz',
     &   'qlib_gra_D16emt_pahi_3.000.gz'/

! D16 graphite, PAH0, a_t=20A, f_gmin=0.01
!               EMT used for graphite
      DATA (FILENAME(23,JBA),JBA=-7,7)/
     &   'qlib_gra_D16emt_pahn_0.333.gz',
     &   'qlib_gra_D16emt_pahn_0.400.gz',
     &   'qlib_gra_D16emt_pahn_0.500.gz',
     &   'qlib_gra_D16emt_pahn_0.556.gz',
     &   'qlib_gra_D16emt_pahn_0.625.gz',
     &   'qlib_gra_D16emt_pahn_0.714.gz',
     &   'qlib_gra_D16emt_pahn_0.833.gz',
     &   'qlib_gra_D16emt_pahn_1.000.gz',
     &   'qlib_gra_D16emt_pahn_1.200.gz',
     &   'qlib_gra_D16emt_pahn_1.400.gz',
     &   'qlib_gra_D16emt_pahn_1.600.gz',
     &   'qlib_gra_D16emt_pahn_1.800.gz',
     &   'qlib_gra_D16emt_pahn_2.000.gz',
     &   'qlib_gra_D16emt_pahn_2.500.gz',
     &   'qlib_gra_D16emt_pahn_3.000.gz'/

! D16 graphite, PAH+, a_t=50A, f_gmin=0.01
!               EMT used for graphite
      DATA (FILENAME(24,JBA),JBA=-7,7)/
     &   'qlib_gra_D16emt_pahib_0.333.gz',
     &   'qlib_gra_D16emt_pahib_0.400.gz',
     &   'qlib_gra_D16emt_pahib_0.500.gz',
     &   'qlib_gra_D16emt_pahib_0.556.gz',
     &   'qlib_gra_D16emt_pahib_0.625.gz',
     &   'qlib_gra_D16emt_pahib_0.714.gz',
     &   'qlib_gra_D16emt_pahib_0.833.gz',
     &   'qlib_gra_D16emt_pahib_1.000.gz',
     &   'qlib_gra_D16emt_pahib_1.200.gz',
     &   'qlib_gra_D16emt_pahib_1.400.gz',
     &   'qlib_gra_D16emt_pahib_1.600.gz',
     &   'qlib_gra_D16emt_pahib_1.800.gz',
     &   'qlib_gra_D16emt_pahib_2.000.gz',
     &   'qlib_gra_D16emt_pahib_2.500.gz',
     &   'qlib_gra_D16emt_pahib_3.000.gz'/

! D16 graphite, PAH0, a_t=50A, f_gmin=0.01
!               EMT used for graphite
      DATA (FILENAME(25,JBA),JBA=-7,7)/
     &   'qlib_gra_D16emt_pahnb_0.333.gz',
     &   'qlib_gra_D16emt_pahnb_0.400.gz',
     &   'qlib_gra_D16emt_pahnb_0.500.gz',
     &   'qlib_gra_D16emt_pahnb_0.556.gz',
     &   'qlib_gra_D16emt_pahnb_0.625.gz',
     &   'qlib_gra_D16emt_pahnb_0.714.gz',
     &   'qlib_gra_D16emt_pahnb_0.833.gz',
     &   'qlib_gra_D16emt_pahnb_1.000.gz',
     &   'qlib_gra_D16emt_pahnb_1.200.gz',
     &   'qlib_gra_D16emt_pahnb_1.400.gz',
     &   'qlib_gra_D16emt_pahnb_1.600.gz',
     &   'qlib_gra_D16emt_pahnb_1.800.gz',
     &   'qlib_gra_D16emt_pahnb_2.000.gz',
     &   'qlib_gra_D16emt_pahnb_2.500.gz',
     &   'qlib_gra_D16emt_pahnb_3.000.gz'/

! D16 graphite, random c axis
      DATA (FILENAME(26,JBA),JBA=-7,7)/
     &   'qlib_gra_D16rc_0.333.gz',
     &   'qlib_gra_D16rc_0.400.gz',
     &   'qlib_gra_D16rc_0.500.gz',
     &   'qlib_gra_D16rc_0.556.gz',
     &   'qlib_gra_D16rc_0.625.gz',
     &   'qlib_gra_D16rc_0.714.gz',
     &   'qlib_gra_D16rc_0.833.gz',
     &   'qlib_gra_D16rc_1.000.gz',
     &   'qlib_gra_D16rc_1.200.gz',
     &   'qlib_gra_D16rc_1.400.gz',
     &   'qlib_gra_D16rc_1.600.gz',
     &   'qlib_gra_D16rc_1.800.gz',
     &   'qlib_gra_D16rc_2.000.gz',
     &   'qlib_gra_D16rc_2.500.gz',
     &   'qlib_gra_D16rc_3.000.gz'/

! D16 graphite, random c axis
      DATA (FILENAME(27,JBA),JBA=-7,7)/
     &   'qlib_gra_D16ca_0.333.gz',
     &   'qlib_gra_D16ca_0.400.gz',
     &   'qlib_gra_D16ca_0.500.gz',
     &   'qlib_gra_D16ca_0.556.gz',
     &   'qlib_gra_D16ca_0.625.gz',
     &   'qlib_gra_D16ca_0.714.gz',
     &   'qlib_gra_D16ca_0.833.gz',
     &   'qlib_gra_D16ca_1.000.gz',
     &   'qlib_gra_D16ca_1.200.gz',
     &   'qlib_gra_D16ca_1.400.gz',
     &   'qlib_gra_D16ca_1.600.gz',
     &   'qlib_gra_D16ca_1.800.gz',
     &   'qlib_gra_D16ca_2.000.gz',
     &   'qlib_gra_D16ca_2.500.gz',
     &   'qlib_gra_D16ca_3.000.gz'/

! D16 graphite, EMT
      DATA (FILENAME(28,JBA),JBA=-7,7)/
     &   'qlib_gra_D16emt_0.333.gz',
     &   'qlib_gra_D16emt_0.400.gz',
     &   'qlib_gra_D16emt_0.500.gz',
     &   'qlib_gra_D16emt_0.556.gz',
     &   'qlib_gra_D16emt_0.625.gz',
     &   'qlib_gra_D16emt_0.714.gz',
     &   'qlib_gra_D16emt_0.833.gz',
     &   'qlib_gra_D16emt_1.000.gz',
     &   'qlib_gra_D16emt_1.200.gz',
     &   'qlib_gra_D16emt_1.400.gz',
     &   'qlib_gra_D16emt_1.600.gz',
     &   'qlib_gra_D16emt_1.800.gz',
     &   'qlib_gra_D16emt_2.000.gz',
     &   'qlib_gra_D16emt_2.500.gz',
     &   'qlib_gra_D16emt_3.000.gz'/

! corresponding .dat files for reorganized output

      DATA (FILENAME2(1,JBA),JBA=-7,7)/
     &  'qlib_silicate_0.333.dat',
     &'  qlib_silicate_0.400.dat',
     &  'qlib_silicate_0.500.dat',
     &  'qlib_silicate_0.556.dat',
     &  'qlib_silicate_0.625.dat',
     &  'qlib_silicate_0.714.dat',
     &  'qlib_silicate_0.833.dat',
     &  'qlib_silicate_1.000.dat',
     &  'qlib_silicate_1.200.dat',
     &  'qlib_silicate_1.400.dat',
     &  'qlib_silicate_1.600.dat',
     &  'qlib_silicate_1.800.dat',
     &  'qlib_silicate_2.000.dat',
     &  'qlib_silicate_2.500.dat',
     &  'qlib_silicate_3.000.dat'/

      DATA (FILENAME2(2,JBA),JBA=-7,7)/
     &   'qlib_gra-pahi_0.333.dat',
     &   'qlib_gra-pahi_0.400.dat',
     &   'qlib_gra-pahi_0.500.dat',
     &   'qlib_gra-pahi_0.556.dat',
     &   'qlib_gra-pahi_0.625.dat',
     &   'qlib_gra-pahi_0.714.dat',
     &   'qlib_gra-pahi_0.833.dat',
     &   'qlib_gra-pahi_1.000.dat',
     &   'qlib_gra-pahi_1.200.dat',
     &   'qlib_gra-pahi_1.400.dat',
     &   'qlib_gra-pahi_1.600.dat',
     &   'qlib_gra-pahi_1.800.dat',
     &   'qlib_gra-pahi_2.000.dat',
     &   'qlib_gra-pahi_2.500.dat',
     &   'qlib_gra-pahi_3.000.dat'/

      DATA (FILENAME2(3,JBA),JBA=-7,7)/
     &   'qlib_gra-pahn_0.333.dat',
     &   'qlib_gra-pahn_0.400.dat',
     &   'qlib_gra-pahn_0.500.dat',
     &   'qlib_gra-pahn_0.556.dat',
     &   'qlib_gra-pahn_0.625.dat',
     &   'qlib_gra-pahn_0.714.dat',
     &   'qlib_gra-pahn_0.833.dat',
     &   'qlib_gra-pahn_1.000.dat',
     &   'qlib_gra-pahn_1.200.dat',
     &   'qlib_gra-pahn_1.400.dat',
     &   'qlib_gra-pahn_1.600.dat',
     &   'qlib_gra-pahn_1.800.dat',
     &   'qlib_gra-pahn_2.000.dat',
     &   'qlib_gra-pahn_2.500.dat',
     &   'qlib_gra-pahn_3.000.dat'/

      DATA (FILENAME2(4,JBA),JBA=-7,7)/
     &   'qlib_sil_HD14_0.333.dat',
     &   'qlib_sil_HD14_0.400.dat',
     &   'qlib_sil_HD14_0.500.dat',
     &   'qlib_sil_HD14_0.556.dat',
     &   'qlib_sil_HD14_0.625.dat',
     &   'qlib_sil_HD14_0.714.dat',
     &   'qlib_sil_HD14_0.833.dat',
     &   'qlib_sil_HD14_1.000.dat',
     &   'qlib_sil_HD14_1.200.dat',
     &   'qlib_sil_HD14_1.400.dat',
     &   'qlib_sil_HD14_1.600.dat',
     &   'qlib_sil_HD14_1.800.dat',
     &   'qlib_sil_HD14_2.000.dat',
     &   'qlib_sil_HD14_2.500.dat',
     &   'qlib_sil_HD14_3.000.dat'/

      DATA (FILENAME2(5,JBA),JBA=-7,7)/
     &   'qlib_sil_HD14_Fe0.02_0.333.dat',
     &   'qlib_sil_HD14_Fe0.02_0.400.dat',
     &   'qlib_sil_HD14_Fe0.02_0.500.dat',
     &   'qlib_sil_HD14_Fe0.02_0.556.dat',
     &   'qlib_sil_HD14_Fe0.02_0.625.dat',
     &   'qlib_sil_HD14_Fe0.02_0.714.dat',
     &   'qlib_sil_HD14_Fe0.02_0.833.dat',
     &   'qlib_sil_HD14_Fe0.02_1.000.dat',
     &   'qlib_sil_HD14_Fe0.02_1.200.dat',
     &   'qlib_sil_HD14_Fe0.02_1.400.dat',
     &   'qlib_sil_HD14_Fe0.02_1.600.dat',
     &   'qlib_sil_HD14_Fe0.02_1.800.dat',
     &   'qlib_sil_HD14_Fe0.02_2.000.dat',
     &   'qlib_sil_HD14_Fe0.02_2.500.dat',
     &   'qlib_sil_HD14_Fe0.02_3.000.dat'/
 
       DATA (FILENAME2(6,JBA),JBA=-7,7)/
     &   'qlib_sil_HD14_Fe0.05_0.333.dat',
     &   'qlib_sil_HD14_Fe0.05_0.400.dat',
     &   'qlib_sil_HD14_Fe0.05_0.500.dat',
     &   'qlib_sil_HD14_Fe0.05_0.556.dat',
     &   'qlib_sil_HD14_Fe0.05_0.625.dat',
     &   'qlib_sil_HD14_Fe0.05_0.714.dat',
     &   'qlib_sil_HD14_Fe0.05_0.833.dat',
     &   'qlib_sil_HD14_Fe0.05_1.000.dat',
     &   'qlib_sil_HD14_Fe0.05_1.200.dat',
     &   'qlib_sil_HD14_Fe0.05_1.400.dat',
     &   'qlib_sil_HD14_Fe0.05_1.600.dat',
     &   'qlib_sil_HD14_Fe0.05_1.800.dat',
     &   'qlib_sil_HD14_Fe0.05_2.000.dat',
     &   'qlib_sil_HD14_Fe0.05_2.500.dat',
     &   'qlib_sil_HD14_Fe0.05_3.000.dat'/

       DATA (FILENAME2(7,JBA),JBA=-7,7)/
     &   'qlib_sil_HD14_Fe0.10_0.333.dat',
     &   'qlib_sil_HD14_Fe0.10_0.400.dat',
     &   'qlib_sil_HD14_Fe0.10_0.500.dat',
     &   'qlib_sil_HD14_Fe0.10_0.556.dat',
     &   'qlib_sil_HD14_Fe0.10_0.625.dat',
     &   'qlib_sil_HD14_Fe0.10_0.714.dat',
     &   'qlib_sil_HD14_Fe0.10_0.833.dat',
     &   'qlib_sil_HD14_Fe0.10_1.000.dat',
     &   'qlib_sil_HD14_Fe0.10_1.200.dat',
     &   'qlib_sil_HD14_Fe0.10_1.400.dat',
     &   'qlib_sil_HD14_Fe0.10_1.600.dat',
     &   'qlib_sil_HD14_Fe0.10_1.800.dat',
     &   'qlib_sil_HD14_Fe0.10_2.000.dat',
     &   'qlib_sil_HD14_Fe0.10_2.500.dat',
     &   'qlib_sil_HD14_Fe0.10_3.000.dat'/

       DATA (FILENAME2(8,JBA),JBA=-7,7)/
     &   'qlib_sil_HD14_Fe0.20_0.333.dat',
     &   'qlib_sil_HD14_Fe0.20_0.400.dat',
     &   'qlib_sil_HD14_Fe0.20_0.500.dat',
     &   'qlib_sil_HD14_Fe0.20_0.556.dat',
     &   'qlib_sil_HD14_Fe0.20_0.625.dat',
     &   'qlib_sil_HD14_Fe0.20_0.714.dat',
     &   'qlib_sil_HD14_Fe0.20_0.833.dat',
     &   'qlib_sil_HD14_Fe0.20_1.000.dat',
     &   'qlib_sil_HD14_Fe0.20_1.200.dat',
     &   'qlib_sil_HD14_Fe0.20_1.400.dat',
     &   'qlib_sil_HD14_Fe0.20_1.600.dat',
     &   'qlib_sil_HD14_Fe0.20_1.800.dat',
     &   'qlib_sil_HD14_Fe0.20_2.000.dat',
     &   'qlib_sil_HD14_Fe0.20_2.500.dat',
     &   'qlib_sil_HD14_Fe0.20_3.000.dat'/

      DATA (FILENAME2(9,JBA),JBA=-7,7)/
     &   'qlib_DH13_Fe_0.333.dat',
     &   'qlib_DH13_Fe_0.400.dat',
     &   'qlib_DH13_Fe_0.500.dat',
     &   'qlib_DH13_Fe_0.556.dat',
     &   'qlib_DH13_Fe_0.625.dat',
     &   'qlib_DH13_Fe_0.714.dat',
     &   'qlib_DH13_Fe_0.833.dat',
     &   'qlib_DH13_Fe_1.000.dat',
     &   'qlib_DH13_Fe_1.200.dat',
     &   'qlib_DH13_Fe_1.400.dat',
     &   'qlib_DH13_Fe_1.600.dat',
     &   'qlib_DH13_Fe_1.800.dat',
     &   'qlib_DH13_Fe_2.000.dat',
     &   'qlib_DH13_Fe_2.500.dat',
     &   'qlib_DH13_Fe_3.000.dat'/

      DATA (FILENAME2(10,JBA),JBA=-7,7)/
     &   'qlib_gra-pahi_HD14_0.333.dat',
     &   'qlib_gra-pahi_HD14_0.400.dat',
     &   'qlib_gra-pahi_HD14_0.500.dat',
     &   'qlib_gra-pahi_HD14_0.556.dat',
     &   'qlib_gra-pahi_HD14_0.625.dat',
     &   'qlib_gra-pahi_HD14_0.714.dat',
     &   'qlib_gra-pahi_HD14_0.833.dat',
     &   'qlib_gra-pahi_HD14_1.000.dat',
     &   'qlib_gra-pahi_HD14_1.200.dat',
     &   'qlib_gra-pahi_HD14_1.400.dat',
     &   'qlib_gra-pahi_HD14_1.600.dat',
     &   'qlib_gra-pahi_HD14_1.800.dat',
     &   'qlib_gra-pahi_HD14_2.000.dat',
     &   'qlib_gra-pahi_HD14_2.500.dat',
     &   'qlib_gra-pahi_HD14_3.000.dat'/

      DATA (FILENAME2(11,JBA),JBA=-7,7)/
     &   'qlib_gra-pahn_HD14_0.333.dat',
     &   'qlib_gra-pahn_HD14_0.400.dat',
     &   'qlib_gra-pahn_HD14_0.500.dat',
     &   'qlib_gra-pahn_HD14_0.556.dat',
     &   'qlib_gra-pahn_HD14_0.625.dat',
     &   'qlib_gra-pahn_HD14_0.714.dat',
     &   'qlib_gra-pahn_HD14_0.833.dat',
     &   'qlib_gra-pahn_HD14_1.000.dat',
     &   'qlib_gra-pahn_HD14_1.200.dat',
     &   'qlib_gra-pahn_HD14_1.400.dat',
     &   'qlib_gra-pahn_HD14_1.600.dat',
     &   'qlib_gra-pahn_HD14_1.800.dat',
     &   'qlib_gra-pahn_HD14_2.000.dat',
     &   'qlib_gra-pahn_HD14_2.500.dat',
     &   'qlib_gra-pahn_HD14_3.000.dat'/

 
      DATA (FILENAME2(12,JBA),JBA=-7,7)/
     &   'qlib_gra-pahi_HD14b_0.333.dat',
     &   'qlib_gra-pahi_HD14b_0.400.dat',
     &   'qlib_gra-pahi_HD14b_0.500.dat',
     &   'qlib_gra-pahi_HD14b_0.556.dat',
     &   'qlib_gra-pahi_HD14b_0.625.dat',
     &   'qlib_gra-pahi_HD14b_0.714.dat',
     &   'qlib_gra-pahi_HD14b_0.833.dat',
     &   'qlib_gra-pahi_HD14b_1.000.dat',
     &   'qlib_gra-pahi_HD14b_1.200.dat',
     &   'qlib_gra-pahi_HD14b_1.400.dat',
     &   'qlib_gra-pahi_HD14b_1.600.dat',
     &   'qlib_gra-pahi_HD14b_1.800.dat',
     &   'qlib_gra-pahi_HD14b_2.000.dat',
     &   'qlib_gra-pahi_HD14b_2.500.dat',
     &   'qlib_gra-pahi_HD14b_3.000.dat'/

      DATA (FILENAME2(13,JBA),JBA=-7,7)/
     &   'qlib_gra-pahn_HD14b_0.333.dat',
     &   'qlib_gra-pahn_HD14b_0.400.dat',
     &   'qlib_gra-pahn_HD14b_0.500.dat',
     &   'qlib_gra-pahn_HD14b_0.556.dat',
     &   'qlib_gra-pahn_HD14b_0.625.dat',
     &   'qlib_gra-pahn_HD14b_0.714.dat',
     &   'qlib_gra-pahn_HD14b_0.833.dat',
     &   'qlib_gra-pahn_HD14b_1.000.dat',
     &   'qlib_gra-pahn_HD14b_1.200.dat',
     &   'qlib_gra-pahn_HD14b_1.400.dat',
     &   'qlib_gra-pahn_HD14b_1.600.dat',
     &   'qlib_gra-pahn_HD14b_1.800.dat',
     &   'qlib_gra-pahn_HD14b_2.000.dat',
     &   'qlib_gra-pahn_HD14b_2.500.dat',
     &   'qlib_gra-pahn_HD14b_3.000.dat'/

! D16 graphite, PAH+, a_t=20A, f_gmin=0.01
!               randomly-oriented c axis
      DATA (FILENAME2(14,JBA),JBA=-7,7)/
     &   'qlib_gra_D16rc_pahi_0.333.dat',
     &   'qlib_gra_D16rc_pahi_0.400.dat',
     &   'qlib_gra_D16rc_pahi_0.500.dat',
     &   'qlib_gra_D16rc_pahi_0.556.dat',
     &   'qlib_gra_D16rc_pahi_0.625.dat',
     &   'qlib_gra_D16rc_pahi_0.714.dat',
     &   'qlib_gra_D16rc_pahi_0.833.dat',
     &   'qlib_gra_D16rc_pahi_1.000.dat',
     &   'qlib_gra_D16rc_pahi_1.200.dat',
     &   'qlib_gra_D16rc_pahi_1.400.dat',
     &   'qlib_gra_D16rc_pahi_1.600.dat',
     &   'qlib_gra_D16rc_pahi_1.800.dat',
     &   'qlib_gra_D16rc_pahi_2.000.dat',
     &   'qlib_gra_D16rc_pahi_2.500.dat',
     &   'qlib_gra_D16rc_pahi_3.000.dat'/

! D16 graphite, PAH0, a_t=20A, f_gmin=0.01
!               randomly-oriented c axis
      DATA (FILENAME2(15,JBA),JBA=-7,7)/
     &   'qlib_gra_D16rc_pahn_0.333.dat',
     &   'qlib_gra_D16rc_pahn_0.400.dat',
     &   'qlib_gra_D16rc_pahn_0.500.dat',
     &   'qlib_gra_D16rc_pahn_0.556.dat',
     &   'qlib_gra_D16rc_pahn_0.625.dat',
     &   'qlib_gra_D16rc_pahn_0.714.dat',
     &   'qlib_gra_D16rc_pahn_0.833.dat',
     &   'qlib_gra_D16rc_pahn_1.000.dat',
     &   'qlib_gra_D16rc_pahn_1.200.dat',
     &   'qlib_gra_D16rc_pahn_1.400.dat',
     &   'qlib_gra_D16rc_pahn_1.600.dat',
     &   'qlib_gra_D16rc_pahn_1.800.dat',
     &   'qlib_gra_D16rc_pahn_2.000.dat',
     &   'qlib_gra_D16rc_pahn_2.500.dat',
     &   'qlib_gra_D16rc_pahn_3.000.dat'/

! D16 graphite, PAH+, a_t=50A, f_gmin=0.01
!               randomly-oriented c axis
      DATA (FILENAME2(16,JBA),JBA=-7,7)/
     &   'qlib_gra_D16rc_pahib_0.333.dat',
     &   'qlib_gra_D16rc_pahib_0.400.dat',
     &   'qlib_gra_D16rc_pahib_0.500.dat',
     &   'qlib_gra_D16rc_pahib_0.556.dat',
     &   'qlib_gra_D16rc_pahib_0.625.dat',
     &   'qlib_gra_D16rc_pahib_0.714.dat',
     &   'qlib_gra_D16rc_pahib_0.833.dat',
     &   'qlib_gra_D16rc_pahib_1.000.dat',
     &   'qlib_gra_D16rc_pahib_1.200.dat',
     &   'qlib_gra_D16rc_pahib_1.400.dat',
     &   'qlib_gra_D16rc_pahib_1.600.dat',
     &   'qlib_gra_D16rc_pahib_1.800.dat',
     &   'qlib_gra_D16rc_pahib_2.000.dat',
     &   'qlib_gra_D16rc_pahib_2.500.dat',
     &   'qlib_gra_D16rc_pahib_3.000.dat'/

! D16 graphite, PAH0, a_t=20A, f_gmin=0.01
!               randomly-oriented c axis
      DATA (FILENAME2(17,JBA),JBA=-7,7)/
     &   'qlib_gra_D16rc_pahnb_0.333.dat',
     &   'qlib_gra_D16rc_pahnb_0.400.dat',
     &   'qlib_gra_D16rc_pahnb_0.500.dat',
     &   'qlib_gra_D16rc_pahnb_0.556.dat',
     &   'qlib_gra_D16rc_pahnb_0.625.dat',
     &   'qlib_gra_D16rc_pahnb_0.714.dat',
     &   'qlib_gra_D16rc_pahnb_0.833.dat',
     &   'qlib_gra_D16rc_pahnb_1.000.dat',
     &   'qlib_gra_D16rc_pahnb_1.200.dat',
     &   'qlib_gra_D16rc_pahnb_1.400.dat',
     &   'qlib_gra_D16rc_pahnb_1.600.dat',
     &   'qlib_gra_D16rc_pahnb_1.800.dat',
     &   'qlib_gra_D16rc_pahnb_2.000.dat',
     &   'qlib_gra_D16rc_pahnb_2.500.dat',
     &   'qlib_gra_D16rc_pahnb_3.000.dat'/

! D16 graphite, PAH+, a_t=20A, f_gmin=0.01
!               c axis || symm axis a
      DATA (FILENAME2(18,JBA),JBA=-7,7)/
     &   'qlib_gra_D16ca_pahi_0.333.dat',
     &   'qlib_gra_D16ca_pahi_0.400.dat',
     &   'qlib_gra_D16ca_pahi_0.500.dat',
     &   'qlib_gra_D16ca_pahi_0.556.dat',
     &   'qlib_gra_D16ca_pahi_0.625.dat',
     &   'qlib_gra_D16ca_pahi_0.714.dat',
     &   'qlib_gra_D16ca_pahi_0.833.dat',
     &   'qlib_gra_D16ca_pahi_1.000.dat',
     &   'qlib_gra_D16ca_pahi_1.200.dat',
     &   'qlib_gra_D16ca_pahi_1.400.dat',
     &   'qlib_gra_D16ca_pahi_1.600.dat',
     &   'qlib_gra_D16ca_pahi_1.800.dat',
     &   'qlib_gra_D16ca_pahi_2.000.dat',
     &   'qlib_gra_D16ca_pahi_2.500.dat',
     &   'qlib_gra_D16ca_pahi_3.000.dat'/

! D16 graphite, PAH0, a_t=20A, f_gmin=0.01
!               c axis || symm axis a
      DATA (FILENAME2(19,JBA),JBA=-7,7)/
     &   'qlib_gra_D16ca_pahn_0.333.dat',
     &   'qlib_gra_D16ca_pahn_0.400.dat',
     &   'qlib_gra_D16ca_pahn_0.500.dat',
     &   'qlib_gra_D16ca_pahn_0.556.dat',
     &   'qlib_gra_D16ca_pahn_0.625.dat',
     &   'qlib_gra_D16ca_pahn_0.714.dat',
     &   'qlib_gra_D16ca_pahn_0.833.dat',
     &   'qlib_gra_D16ca_pahn_1.000.dat',
     &   'qlib_gra_D16ca_pahn_1.200.dat',
     &   'qlib_gra_D16ca_pahn_1.400.dat',
     &   'qlib_gra_D16ca_pahn_1.600.dat',
     &   'qlib_gra_D16ca_pahn_1.800.dat',
     &   'qlib_gra_D16ca_pahn_2.000.dat',
     &   'qlib_gra_D16ca_pahn_2.500.dat',
     &   'qlib_gra_D16ca_pahn_3.000.dat'/

! D16 graphite, PAH+, a_t=50A, f_gmin=0.01
!               c axis || symm axis a
      DATA (FILENAME2(20,JBA),JBA=-7,7)/
     &   'qlib_gra_D16ca_pahib_0.333.dat',
     &   'qlib_gra_D16ca_pahib_0.400.dat',
     &   'qlib_gra_D16ca_pahib_0.500.dat',
     &   'qlib_gra_D16ca_pahib_0.556.dat',
     &   'qlib_gra_D16ca_pahib_0.625.dat',
     &   'qlib_gra_D16ca_pahib_0.714.dat',
     &   'qlib_gra_D16ca_pahib_0.833.dat',
     &   'qlib_gra_D16ca_pahib_1.000.dat',
     &   'qlib_gra_D16ca_pahib_1.200.dat',
     &   'qlib_gra_D16ca_pahib_1.400.dat',
     &   'qlib_gra_D16ca_pahib_1.600.dat',
     &   'qlib_gra_D16ca_pahib_1.800.dat',
     &   'qlib_gra_D16ca_pahib_2.000.dat',
     &   'qlib_gra_D16ca_pahib_2.500.dat',
     &   'qlib_gra_D16ca_pahib_3.000.dat'/

! D16 graphite, PAH0, a_t=50A, f_gmin=0.01
!               c axis || symm axis a
      DATA (FILENAME2(21,JBA),JBA=-7,7)/
     &   'qlib_gra_D16ca_pahnb_0.333.dat',
     &   'qlib_gra_D16ca_pahnb_0.400.dat',
     &   'qlib_gra_D16ca_pahnb_0.500.dat',
     &   'qlib_gra_D16ca_pahnb_0.556.dat',
     &   'qlib_gra_D16ca_pahnb_0.625.dat',
     &   'qlib_gra_D16ca_pahnb_0.714.dat',
     &   'qlib_gra_D16ca_pahnb_0.833.dat',
     &   'qlib_gra_D16ca_pahnb_1.000.dat',
     &   'qlib_gra_D16ca_pahnb_1.200.dat',
     &   'qlib_gra_D16ca_pahnb_1.400.dat',
     &   'qlib_gra_D16ca_pahnb_1.600.dat',
     &   'qlib_gra_D16ca_pahnb_1.800.dat',
     &   'qlib_gra_D16ca_pahnb_2.000.dat',
     &   'qlib_gra_D16ca_pahnb_2.500.dat',
     &   'qlib_gra_D16ca_pahnb_3.000.dat'/

! D16 graphite, PAH+, a_t=20A, f_gmin=0.01
!               EMT used for graphite
      DATA (FILENAME2(22,JBA),JBA=-7,7)/
     &   'qlib_gra_D16emt_pahi_0.333.dat',
     &   'qlib_gra_D16emt_pahi_0.400.dat',
     &   'qlib_gra_D16emt_pahi_0.500.dat',
     &   'qlib_gra_D16emt_pahi_0.556.dat',
     &   'qlib_gra_D16emt_pahi_0.625.dat',
     &   'qlib_gra_D16emt_pahi_0.714.dat',
     &   'qlib_gra_D16emt_pahi_0.833.dat',
     &   'qlib_gra_D16emt_pahi_1.000.dat',
     &   'qlib_gra_D16emt_pahi_1.200.dat',
     &   'qlib_gra_D16emt_pahi_1.400.dat',
     &   'qlib_gra_D16emt_pahi_1.600.dat',
     &   'qlib_gra_D16emt_pahi_1.800.dat',
     &   'qlib_gra_D16emt_pahi_2.000.dat',
     &   'qlib_gra_D16emt_pahi_2.500.dat',
     &   'qlib_gra_D16emt_pahi_3.000.dat'/

! D16 graphite, PAH0, a_t=20A, f_gmin=0.01
!               EMT used for graphite
      DATA (FILENAME2(23,JBA),JBA=-7,7)/
     &   'qlib_gra_D16emt_pahn_0.333.dat',
     &   'qlib_gra_D16emt_pahn_0.400.dat',
     &   'qlib_gra_D16emt_pahn_0.500.dat',
     &   'qlib_gra_D16emt_pahn_0.556.dat',
     &   'qlib_gra_D16emt_pahn_0.625.dat',
     &   'qlib_gra_D16emt_pahn_0.714.dat',
     &   'qlib_gra_D16emt_pahn_0.833.dat',
     &   'qlib_gra_D16emt_pahn_1.000.dat',
     &   'qlib_gra_D16emt_pahn_1.200.dat',
     &   'qlib_gra_D16emt_pahn_1.400.dat',
     &   'qlib_gra_D16emt_pahn_1.600.dat',
     &   'qlib_gra_D16emt_pahn_1.800.dat',
     &   'qlib_gra_D16emt_pahn_2.000.dat',
     &   'qlib_gra_D16emt_pahn_2.500.dat',
     &   'qlib_gra_D16emt_pahn_3.000.dat'/

! D16 graphite, PAH+, a_t=50A, f_gmin=0.01
!               EMT used for graphite
      DATA (FILENAME2(24,JBA),JBA=-7,7)/
     &   'qlib_gra_D16emt_pahib_0.333.dat',
     &   'qlib_gra_D16emt_pahib_0.400.dat',
     &   'qlib_gra_D16emt_pahib_0.500.dat',
     &   'qlib_gra_D16emt_pahib_0.556.dat',
     &   'qlib_gra_D16emt_pahib_0.625.dat',
     &   'qlib_gra_D16emt_pahib_0.714.dat',
     &   'qlib_gra_D16emt_pahib_0.833.dat',
     &   'qlib_gra_D16emt_pahib_1.000.dat',
     &   'qlib_gra_D16emt_pahib_1.200.dat',
     &   'qlib_gra_D16emt_pahib_1.400.dat',
     &   'qlib_gra_D16emt_pahib_1.600.dat',
     &   'qlib_gra_D16emt_pahib_1.800.dat',
     &   'qlib_gra_D16emt_pahib_2.000.dat',
     &   'qlib_gra_D16emt_pahib_2.500.dat',
     &   'qlib_gra_D16emt_pahib_3.000.dat'/

! D16 graphite, PAH0, a_t=50A, f_gmin=0.01
!               EMT used for graphite
      DATA (FILENAME2(25,JBA),JBA=-7,7)/
     &   'qlib_gra_D16emt_pahnb_0.333.dat',
     &   'qlib_gra_D16emt_pahnb_0.400.dat',
     &   'qlib_gra_D16emt_pahnb_0.500.dat',
     &   'qlib_gra_D16emt_pahnb_0.556.dat',
     &   'qlib_gra_D16emt_pahnb_0.625.dat',
     &   'qlib_gra_D16emt_pahnb_0.714.dat',
     &   'qlib_gra_D16emt_pahnb_0.833.dat',
     &   'qlib_gra_D16emt_pahnb_1.000.dat',
     &   'qlib_gra_D16emt_pahnb_1.200.dat',
     &   'qlib_gra_D16emt_pahnb_1.400.dat',
     &   'qlib_gra_D16emt_pahnb_1.600.dat',
     &   'qlib_gra_D16emt_pahnb_1.800.dat',
     &   'qlib_gra_D16emt_pahnb_2.000.dat',
     &   'qlib_gra_D16emt_pahnb_2.500.dat',
     &   'qlib_gra_D16emt_pahnb_3.000.dat'/

! D16 graphite, random c axis
      DATA (FILENAME2(26,JBA),JBA=-7,7)/
     &   'qlib_gra_D16rc_0.333.dat',
     &   'qlib_gra_D16rc_0.400.dat',
     &   'qlib_gra_D16rc_0.500.dat',
     &   'qlib_gra_D16rc_0.556.dat',
     &   'qlib_gra_D16rc_0.625.dat',
     &   'qlib_gra_D16rc_0.714.dat',
     &   'qlib_gra_D16rc_0.833.dat',
     &   'qlib_gra_D16rc_1.000.dat',
     &   'qlib_gra_D16rc_1.200.dat',
     &   'qlib_gra_D16rc_1.400.dat',
     &   'qlib_gra_D16rc_1.600.dat',
     &   'qlib_gra_D16rc_1.800.dat',
     &   'qlib_gra_D16rc_2.000.dat',
     &   'qlib_gra_D16rc_2.500.dat',
     &   'qlib_gra_D16rc_3.000.dat'/

! D16 graphite, random c axis
      DATA (FILENAME2(27,JBA),JBA=-7,7)/
     &   'qlib_gra_D16ca_0.333.dat',
     &   'qlib_gra_D16ca_0.400.dat',
     &   'qlib_gra_D16ca_0.500.dat',
     &   'qlib_gra_D16ca_0.556.dat',
     &   'qlib_gra_D16ca_0.625.dat',
     &   'qlib_gra_D16ca_0.714.dat',
     &   'qlib_gra_D16ca_0.833.dat',
     &   'qlib_gra_D16ca_1.000.dat',
     &   'qlib_gra_D16ca_1.200.dat',
     &   'qlib_gra_D16ca_1.400.dat',
     &   'qlib_gra_D16ca_1.600.dat',
     &   'qlib_gra_D16ca_1.800.dat',
     &   'qlib_gra_D16ca_2.000.dat',
     &   'qlib_gra_D16ca_2.500.dat',
     &   'qlib_gra_D16ca_3.000.dat'/

! D16 graphite, EMT
      DATA (FILENAME2(28,JBA),JBA=-7,7)/
     &   'qlib_gra_D16emt_0.333.dat',
     &   'qlib_gra_D16emt_0.400.dat',
     &   'qlib_gra_D16emt_0.500.dat',
     &   'qlib_gra_D16emt_0.556.dat',
     &   'qlib_gra_D16emt_0.625.dat',
     &   'qlib_gra_D16emt_0.714.dat',
     &   'qlib_gra_D16emt_0.833.dat',
     &   'qlib_gra_D16emt_1.000.dat',
     &   'qlib_gra_D16emt_1.200.dat',
     &   'qlib_gra_D16emt_1.400.dat',
     &   'qlib_gra_D16emt_1.600.dat',
     &   'qlib_gra_D16emt_1.800.dat',
     &   'qlib_gra_D16emt_2.000.dat',
     &   'qlib_gra_D16emt_2.500.dat',
     &   'qlib_gra_D16emt_3.000.dat'/

! DH21 Astrodust

      DATA (CHPO(JPO),JPO=0,9)/
     &   'qlib_DH21Ad_P0.00','qlib_DH21Ad_P0.10','qlib_DH21Ad_P0.20',
     &   'qlib_DH21Ad_P0.30','qlib_DH21Ad_P0.40','qlib_DH21Ad_P0.50',
     &   'qlib_DH21Ad_P0.60','qlib_DH21Ad_P0.70','qlib_DH21Ad_P0.80',
     &   'qlib_DH21Ad_P0.90'/
      DATA (CHPO2(JPO),JPO=0,9)/
     &   'q_DH21Ad_P0.00','q_DH21Ad_P0.10','q_DH21Ad_P0.20',
     &   'q_DH21Ad_P0.30','q_DH21Ad_P0.40','q_DH21Ad_P0.50',
     &   'q_DH21Ad_P0.60','q_DH21Ad_P0.70','q_DH21Ad_P0.80',
     &   'q_DH21Ad_P0.90'/
      DATA (CHFE(JFE),JFE=0,5)/
     &   '_Fe0.00','_Fe0.10','_Fe0.20','_Fe0.30','_Fe0.40','_Fe0.50'/
      DATA (CHBA(JBA),JBA=-7,7)/
     &   '_0.333','_0.400','_0.500','_0.556','_0.625','_0.714','_0.833',
     &   '_1.000','_1.200','_1.400','_1.600','_1.800','_2.000','_2.500',
     &   '_3.000'/

! JC from 29 to 28 + 10*6 = 88

      DO JPO=0,9
         DO JFE=0,5
            IDX=JFE+JPO*6
            JC=29+IDX
            DO JBA=-7,7
               FILENAME(JC,JBA)=
     &            CHPO(JPO)//CHFE(JFE)//CHBA(JBA)//'.gz'

! filename2 is for output

               FILENAME2(JC,JBA)=
     &            CHPO2(JPO)//CHFE(JFE)//CHBA(JBA)//'.dat'
            ENDDO
         ENDDO
      ENDDO

! JC from 89 to 91

      DO JBA=-7,7
         FILENAME(89,JBA)='qlib_gra_D16emt2'//CHBA(JBA)//'.gz'
         FILENAME(90,JBA)='qlib_gra_D16emt2_pahib'//CHBA(JBA)//'.gz'
         FILENAME(91,JBA)='qlib_gra_D16emt2_pahnb'//CHBA(JBA)//'.gz'

! 2021.02.01 (BTD) fix: add following 3 lines

         FILENAME2(89,JBA)='qlib_gra_D16emt2'//CHBA(JBA)//'.dat'
         FILENAME2(90,JBA)='qlib_gra_D16emt2_pahib'//CHBA(JBA)//'.dat'
         FILENAME2(91,JBA)='qlib_gra_D16emt2_pahnb'//CHBA(JBA)//'.dat'

      ENDDO

!*** diagnostic
!      write(0,fmt='(a,i4,a,1pe10.3,a,1pe10.3)')
!     &   'readqlib_v19 ckpt 0, icase=',icase,
!     &   ' bovera=',bovera,' wave=',wave
!***
      IF(NCASE.EQ.0)THEN
         PATH='../data/'

         DO JC=1,NCTMAX
            DO JS=-7,7
               ICT(JC,JS)=0
            ENDDO
         ENDDO
      ENDIF

! given BOVERA, determine ISHAP

      IF(ABS(BOVERA-0.333).LT.0.001)THEN
         ISHAP=-7
      ELSEIF(ABS(BOVERA-0.400).LT.0.001)THEN
         ISHAP=-6
      ELSEIF(ABS(BOVERA-0.500).LT.0.001)THEN
         ISHAP=-5
      ELSEIF(ABS(BOVERA-0.556).LT.0.001)THEN
         ISHAP=-4
      ELSEIF(ABS(BOVERA-0.625).LT.0.001)THEN
         ISHAP=-3
      ELSEIF(ABS(BOVERA-0.714).LT.0.001)THEN
         ISHAP=-2
      ELSEIF(ABS(BOVERA-0.833).LT.0.001)THEN
         ISHAP=-1
      ELSEIF(ABS(BOVERA-1.000).LT.0.001)THEN
         ISHAP=0
      ELSEIF(ABS(BOVERA-1.200).LT.0.001)THEN
         ISHAP=1
      ELSEIF(ABS(BOVERA-1.400).LT.0.001)THEN
         ISHAP=2
      ELSEIF(ABS(BOVERA-1.600).LT.0.001)THEN
         ISHAP=3
      ELSEIF(ABS(BOVERA-1.800).LT.0.001)THEN
         ISHAP=4
      ELSEIF(ABS(BOVERA-2.000).LT.0.001)THEN
         ISHAP=5
      ELSEIF(ABS(BOVERA-2.500).LT.0.001)THEN
         ISHAP=6
      ELSEIF(ABS(BOVERA-3.000).LT.0.001)THEN
         ISHAP=7
      ELSE
         WRITE(0,FMT='(A,F10.5)')
     &      'readqlib_v19 ckpt 1 fatal error: invalid BOVERA=',BOVERA
         STOP
      ENDIF

! initialization if necessary:

      IF(ICT(ICASE,ISHAP).EQ.0)THEN
         NCASE=NCASE+1

! verify that memory is available to store this new case

         IF(NCASE.GT.NCMAX)THEN

! free up memory for new case
! discard data now stored for JC=1

!*** diagnostic
!            WRITE(0,FMT='(A,A,I2,A,A,I2)')'readqlib_v19 ckpt 2:',
!     &         ' called for ',NCASE,' different compositions',
!     &           ' but memory allocated only for NCMAX=',NCMAX
!            write(0,fmt='(a,i10)')'icase0(1)=',icase0(1)
!            write(0,fmt='(a,i10)')'ishap0(1)=',ishap0(1)
!***
            ICT(ICASE0(1),ISHAP0(1))=0
            DO JC=2,NCMAX
!*** diagnostic
!               WRITE(0,FMT='(A,I3)')'readqlib_v19 ckpt 3, JC=',JC
!***
               NRAD(JC-1)=NRAD(JC)
               NWAV(JC-1)=NWAV(JC)
               AXRAT(JC-1)=AXRAT(JC)
               DO JR=0,NRAD(JC)
                  RAD(JR,JC-1)=RAD(JR,JC)
               ENDDO
               DO JW=0,NWAV(JC)
                  WAV(JW,JC-1)=WAV(JW,JC)
               ENDDO
!*** diagnostic
!               WRITE(0,FMT='(A,I3)')'readqlib_v19 ckpt 4, JC=',JC
!***
               DO JW=0,NWAV(JC)
                  DO JR=0,NRAD(JC)
                     DO JORI=1,3
                        QTABS(JORI,JR,JW,JC-1)=QTABS(JORI,JR,JW,JC)
                        QTEXT(JORI,JR,JW,JC-1)=QTEXT(JORI,JR,JW,JC)
                        QTSCA(JORI,JR,JW,JC-1)=QTSCA(JORI,JR,JW,JC)
                     ENDDO
                  ENDDO
               ENDDO
               ICASE0(JC-1)=ICASE0(JC)
               ISHAP0(JC-1)=ISHAP0(JC)
               ICT(ICASE0(JC-1),ISHAP0(JC-1))=JC-1
            ENDDO
!*** diagnostic
!           WRITE(0,FMT='(A)')'readqlib_v19 ckpt 5'
!***
            NCASE=NCMAX
         ENDIF
         ICT(ICASE,ISHAP)=NCASE
         ICASE0(NCASE)=ICASE
         ISHAP0(NCASE)=ISHAP

!*** diagnostic
!         WRITE(0,FMT='(A,A)')'readqlib_v19 ckpt 6: about to read file=',
!     &      FILENAME(ICASE,ISHAP)
!***
         
! pass filename to COMMON
         FLNAME=FILENAME(ICASE,ISHAP)

         CALL SYSTEM_CLOCK(COUNT=C1,COUNT_RATE=CRATE,COUNT_MAX=CMAX)
!         WRITE(0,*)'COUNT=',C1
!         WRITE(0,*)'COUNT_RATE=',CRATE
!         WRITE(0,*)'COUNT_MAX=',CMAX
         J1=MOD(C1,10)
         J2=MOD(INT(C1/10),10)
         J3=MOD(INT(C1/100),10)
         J4=MOD(INT(C1/1000),10)
         J5=MOD(INT(C1/10000),10)
         J6=MOD(INT(C1/100000),10)
         J7=MOD(INT(C1/1000000),10)
         J8=MOD(INT(C1/10000000),10)
         TMPFILE='tmpfile'//C(J1)//C(J2)//C(J3)//C(J4)//C(J5)//C(J6)//
     &      C(J7)//C(J8)
         CALL SYSTEM('gunzip -c '
     &               //PATH//FILENAME(ICASE,ISHAP)//' > '//TMPFILE)
         OPEN(UNIT=3,FILE=TMPFILE,FORM='FORMATTED',
     &      STATUS='OLD',ERR=9000)
         WRITE(0,FMT='(A,A)')'readqlib_v19 ckpt 7: opened file= ',
     &        PATH//FILENAME(ICASE,ISHAP)
         WRITE(0,FMT='(A,A)')'                using TMPFILE=',TMPFILE
!--- skip one line (descriptor)
         READ(3,*)
!------------------------------
         READ(3,*)NRAD(NCASE),NWAV(NCASE)

!*** diagnostic
         write(0,fmt='(a,i6)')'readqlib_v19 ckpt 8 NRAD=',NRAD(NCASE)
         write(0,fmt='(a,i6)')'                    NWAV=',NWAV(NCASE)
!***

! verify that stored table is compatible with memory allocation

         IF(NRAD(NCASE).GT.NRMAX)THEN
            WRITE(0,FMT='(A,A,A,I6,A,I6)')'Fatal error in READQLIB:',
     &         FILENAME,' has NRAD=',NRAD(NCASE),
     &         ' but READQLIB has only allocated for NRAD=',NRMAX
            STOP
         ENDIF
         IF(NWAV(NCASE).GT.NWMAX)THEN
            WRITE(0,FMT='(A,A,A,I6,A,I6)')'Fatal error in READQLIB:',
     &         FILENAME,' has NWAV=',NWAV(NCASE),
     &         ' but READQLIB has only allocated for NWAV=',NWMAX
            STOP
         ENDIF

!*** diagnostic
!         write(0,fmt='(a,i4)')'readqlib_v19 ckpt 9,  ncase=',ncase
!***

! OK: proceed to read

         READ(3,9200)AXRAT(NCASE)
!-- skip one line --
         READ(3,*)
!-------------------
         READ(3,9200)(RAD(JR,NCASE),JR=0,NRAD(NCASE))
!-- skip one line --
         READ(3,*)
!-------------------
         READ(3,9200)(WAV(JW,NCASE),JW=0,NWAV(NCASE))
!*** diagnostic
!         write(0,fmt='(a)')'readqlib_v19 ckpt 10'
!***
!-- skip one line --
         READ(3,*)
!-------------------
         READ(3,9200)(((QTABS(JORI,JR,JW,NCASE),JORI=1,3),
     &                  JR=0,NRAD(NCASE)),JW=0,NWAV(NCASE))
!*** diagnostic
!         write(0,fmt='(a)')'readqlib_v19 ckpt 11'
!***
!-- skip one line --
         READ(3,*)
!-------------------
         READ(3,9200)(((QTEXT(JORI,JR,JW,NCASE),JORI=1,3),
     &                  JR=0,NRAD(NCASE)),JW=0,NWAV(NCASE))
!*** diagnostic
!         write(0,fmt='(a)')'readqlib_v19 ckpt 12'
!***
!-- skip one line --
         READ(3,*)
!-------------------
         READ(3,9200)(((QTSCA(JORI,JR,JW,NCASE),JORI=1,3),
     &                  JR=0,NRAD(NCASE)),JW=0,NWAV(NCASE))
!*** diagnostic
!         write(0,fmt='(a)')'readqlib_v19 ckpt 13'
!***
         CLOSE(3)
!*** diagnostic
!         write(0,fmt='(a,a)')'readqlib_v19 ckpt 14: TMPFILE=',TMPFILE
!***
         CALL SYSTEM('rm '//TMPFILE)
!*** diagnostic
!         write(0,fmt='(a,i4)')'readqlib_v19 ckpt 15, NCASE=',NCASE
!***

         OPEN(UNIT=3,FILE='DH21_aeff',FORM='FORMATTED',STATUS='UNKNOWN',
     &        ERR=9002)
         WRITE(3,FMT='(A)')'DH21Ad...'
         WRITE(3,FMT='(A,A)')'List of radii(um) aeff(jr), jr=0-168,',
     &      ' [40 per decade]'
         WRITE(3,FMT='(1PE9.3,1P168E10.3)')
     &      (RAD(JR,NCASE),JR=0,NRAD(NCASE))
         CLOSE(3)

         OPEN(UNIT=3,FILE='DH21_wave',FORM='FORMATTED',STATUS='UNKNOWN',
     &        ERR=9002)
         WRITE(3,FMT='(A)')'DH21Ad...'
         WRITE(3,FMT='(A,A)')'List of wavelengths(um) wave(jw),',
     &     ' jw=0-1128  [200 per decade]'
         WRITE(3,FMT='(1PE9.3,1P1128E10.3)')
     &      (WAV(JW,NCASE),JW=0,NWAV(NCASE))
         CLOSE(3)

!*** diagnostic
!         write(0,fmt='(a,2i4)')'readqlib_v19 ckpt 16, ICASE,ISHAP=',
!     &      ICASE,ISHAP
!***
         IDX=ICASE-29
         JPO=IDX/6
         JFE=IDX-6*JPO

!*** diagnostic
         write(0,fmt='(a)')'readqlib_v19 ckpt 17'
         write(0,fmt='(a,i6,a,i6)')'   icase=',icase,' ishap=',ishap
         write(0,fmt='(a,a)')
     &      '   about to open file=filename2(icase,ishap)=',
     &      filename2(icase,ishap)
!***

         OPEN(UNIT=3,FILE=FILENAME2(ICASE,ISHAP),FORM='FORMATTED',
     &        STATUS='UNKNOWN',ERR=9002)

!*** diagnostic
!         write(0,fmt='(a)')'readqlib_v19 ckpt 18'
!***

! 6 comment lines:

         WRITE(3,FMT='(2A)')
     &      'Optical cross sections for Draine & Hensley (2021) ',
     &      'Astrodust spheroid with' 
         WRITE(3,FMT='(F5.2,A)')
     &      (0.1*JPO),' = Porosity'
         WRITE(3,FMT='(F5.2,A)')
     &      (0.1*JFE),' = f_Fe = fraction of Fe in metallic form'
         WRITE(3,FMT='(F5.3,2A)')
     &      BOVERA,' = b/a for spheroid ',
     &      '(< 1 for prolate, > 1 for oblate)'
         WRITE(3,FMT='(4A)')
     &      'Q_ext=C_ext/(pi a_eff^2), ',
     &      'Q_abs=C_abs/(pi a_eff^2), ',
     &      'Q_sca=C_sca/(pi a_eff^2) ',
     &      'calculated for'
         WRITE(3,FMT='(5A)')
     &      '3 orientations: ',
     &      'jori=1: k || a;;  ',
     &      'jori=2: k perp a, E || a;  ',
     &      'jori=3: k perp a, E perp a ',
     &      ' [a=spheroid symmetry axis]'
         WRITE(3,FMT='(A,1PE10.3,A,1PE10.3,A)')
     &      ' 169 sizes: jr=0,168 from a_eff =',RAD(0,NCASE),
     &      ' to',RAD(NRAD(NCASE),NCASE),
     &      ' um (see DH21_aeff)'
         WRITE(3,FMT='(A,1PE10.3,A,1PE10.3,A)')
     &      '1129 wavelengths: jw=0,1128 from lambda=',WAV(0,NCASE),
     &      ' to',WAV(NWAV(NCASE),NCASE),
     &      ' um (see DH21_aeff)'
         WRITE(3,FMT='(A)')
     &      'for each wavelength jw, list '

! first plan:
!         WRITE(3,FMT='(A)')
!     &      ' (Qext(jori,jr,jw),jori=1,3),jr=0,168)'
!         WRITE(3,FMT='(A)')
!     &      ' (Qabs(jori,jr,jw),jori=1,3),jr=0,168)'
!         WRITE(3,FMT='(A)')
!     &      ' (Qsca(jori,jr,jw),jori=1,3),jr=0,168)'

! 169 sizes: 3*169-1=506

!         DO JW=0,NWAV(NCASE)
!            WRITE(3,FMT='(1PE9.3,1P506E10.3)')
!     &         ((QTEXT(JORI,JR,JW,NCASE),JR=0,NRAD(NCASE)),JORI=1,3)
!            WRITE(3,FMT='(1PE9.3,1P506E10.3)')
!     &         ((QTABS(JORI,JR,JW,NCASE),JR=0,NRAD(NCASE)),JORI=1,3)
!            WRITE(3,FMT='(1PE9.3,1P506E10.3)')
!     &         ((QTSCA(JORI,JR,JW,NCASE),JR=0,NRAD(NCASE)),JORI=1,3)
!         ENDDO

! second plan:

         WRITE(3,FMT='(A,I4,A,I3,A)')
     &      ' ((Qext(jw,jr,jori),jw=0,',NWAV(NCASE),',jr=0,',
     &      NRAD(NCASE),'),jori=1,3'
         WRITE(3,FMT='(A,I4,A,I3,A)')
     &      ' ((Qabs(jw,jr,jori),jw=0,',NWAV(NCASE),',jr=0,',
     &      NRAD(NCASE),'),jori=1,3'
         WRITE(3,FMT='(A,I4,A,I3,A)')
     &      ' ((Qsca(jw,jr,jori),jw=0,',NWAV(NCASE),',jr=0,',
     &      NRAD(NCASE),'),jori=1,3'

! 169 sizes: 3*169-1=506

         WRITE(3,FMT='(1PE9.3,1P168E10.3)')
     &      (((QTEXT(JORI,JR,JW,NCASE),
     &      JR=0,NRAD(NCASE)),JW=0,NWAV(NCASE)),JORI=1,3)
         WRITE(3,FMT='(1PE9.3,1P168E10.3)')
     &      (((QTABS(JORI,JR,JW,NCASE),
     &      JR=0,NRAD(NCASE)),JW=0,NWAV(NCASE)),JORI=1,3)
         WRITE(3,FMT='(1PE9.3,1P168E10.3)')
     &      (((QTSCA(JORI,JR,JW,NCASE),
     &      JR=0,NRAD(NCASE)),JW=0,NWAV(NCASE)),JORI=1,3)
         CLOSE(3)
!*** diagnostic
!         write(0,fmt='(a,a)')'readqlib_v19 ckpt 19, gzip ',
!     &      FILENAME2(ICASE,ISHAP)
!***
         CALL SYSTEM('gzip -f '//FILENAME2(ICASE,ISHAP))

!*** sanity check
         DO JW=0,NWAV(NCASE)
            DO JR=0,NRAD(NCASE)
               DO JORI=1,3
                  IF(.NOT.(QTABS(JORI,JR,JW,NCASE).GT.0.))THEN
                     WRITE(0,FMT='(A,I2,A,I4,A,I4,A,I2)')
     &               'readqlib_v19 fatal error for jori=',jori,' jr=',
     &               jr,' JW=',JW,' NCASE=',NCASE
                     WRITE(0,*)'QTABS=',QTABS(JORI,JR,JW,NCASE)
                     STOP
                  ENDIF
                  IF(.NOT.(QTEXT(JORI,JR,JW,NCASE).GT.0.))THEN
                     WRITE(0,FMT='(A,I2,A,I4,A,I4,A,I2)')
     &               'readqlib_v19 fatal error for jori=',jori,' jr=',
     &               jr,' JW=',JW,' NCASE=',NCASE
                     WRITE(0,*)'QTEXT=',QTEXT(JORI,JR,JW,NCASE)
                     STOP
                  ENDIF
                  IF(.NOT.(QTSCA(JORI,JR,JW,NCASE).GT.0.))THEN
                     WRITE(0,FMT='(A,I2,A,I4,A,I4,A,I2)')
     &               'readqlib_v19 fatal error for jori=',jori,' jr=',
     &               jr,' JW=',JW,' NCASE=',NCASE
                     WRITE(0,*)'QTSCA=',QTSCA(JORI,JR,JW,NCASE)
                     STOP
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
!*** completed sanity check

!---- completed initialization for new composition ----------------------

      ENDIF

! normal entry point:

!*** diagnostic
!      write(0,*)'readqlib_v19 ckpt 20, wave=',wave
!***

      IC=ICT(ICASE,ISHAP)
      NR=NRAD(IC)
      NW=NWAV(IC)

! data are logarithmically spaced in radius and wavelength
! stored in order of increasing radius and increasing wavelength

      XR=NR*LOG(AEFF/RAD(0,IC))/LOG(RAD(NR,IC)/RAD(0,IC))
      XW=NW*LOG(WAVE/WAV(0,IC))/LOG(WAV(NW,IC)/WAV(0,IC))

!---------------
! issue warning if called for size or wavelength significantly 
! outside range of table
! if desired, could halt calculation at this point to avoid error...
      IF(XR.LT.-0.1.OR.XR.GT.REAL(NR)+0.1)THEN
         WRITE(0,FMT='(A,1PE10.3,A,1PE10.3,A,1PE10.3)')
     &      'WARNING: aeff=',AEFF,' outside allowed range',RAD(0,IC),
     &      ' - ',RAD(NR,IC)
      ENDIF
      IF(XW.LT.-0.1.OR.XW.GT.REAL(NW)+0.1)THEN
         WRITE(0,FMT='(A,1PE10.3,A,1PE10.3,A,1PE10.3)')
     &      'WARNING: wave=',WAVE,' outside allowed range',WAV(0,IC),
     &      ' - ',WAV(NW,IC)
      ENDIF
!----------------
      JR=INT(XR)
      JW=INT(XW)
      IF(JR.LT.0)JR=0
      IF(JR.GT.NR-1)JR=NR-1
      IF(JW.LT.0)JW=0
      IF(JW.GT.NW-1)JW=NW-1
      XR=XR-JR
      XW=XW-JW

! simple weighting factors for interpolation
! interpolate on logarithm to capture power-law behavior

      W00=(1.-XR)*(1.-XW)
      W01=(1.-XR)*XW
      W10=XR*(1.-XW)
      W11=XR*XW

!*** diagnostic
!      write(0,*)'readqlib_v19 ckpt 21, wave=',wave
!      write(0,*)'     jr=',jr
!      write(0,*)'     jw=',jw
!      write(0,*)'     ic=',ic
!***

      DO JORI=1,3
         QABS(JORI)=EXP(W00*LOG(QTABS(JORI,JR,JW,IC))+
     &                  W01*LOG(QTABS(JORI,JR,JW+1,IC))+
     &                  W10*LOG(QTABS(JORI,JR+1,JW,IC))+
     &                  W11*LOG(QTABS(JORI,JR+1,JW+1,IC)))
         QEXT(JORI)=EXP(W00*LOG(QTEXT(JORI,JR,JW,IC))+
     &                  W01*LOG(QTEXT(JORI,JR,JW+1,IC))+
     &                  W10*LOG(QTEXT(JORI,JR+1,JW,IC))+
     &                  W11*LOG(QTEXT(JORI,JR+1,JW+1,IC)))
         QSCA(JORI)=EXP(W00*LOG(QTSCA(JORI,JR,JW,IC))+
     &                  W01*LOG(QTSCA(JORI,JR,JW+1,IC))+
     &                  W10*LOG(QTSCA(JORI,JR+1,JW,IC))+
     &                  W11*LOG(QTSCA(JORI,JR+1,JW+1,IC)))
      ENDDO
!*** diagnostic
!      write(0,*)'readqlib_v19 ckpt 22, wave=',wave
!***
      IF(RQTEST.GT.0)THEN
         WRITE(0,FMT='(A,I3)')'ICASE=',ICASE
         WRITE(0,FMT='(A,1PE10.3)')'aeff=',AEFF
         WRITE(0,FMT='(A,1PE10.3)')'wave=',WAVE
         DO JORI=1,3
            WRITE(0,FMT='(A,I3,A,I3)')'JORI=',JORI,' IC=',IC
            WRITE(0,FMT='(A,1PE10.3)')'QTABS(JORI,JR  ,JW  ,IC)=',
     &                                QTABS(JORI,JR,JW,IC)
            WRITE(0,FMT='(A,1PE10.3)')'QTABS(JORI,JR+1 ,JW  ,IC)=',
     &                                QTABS(JORI,JR+1,JW,IC)
            WRITE(0,FMT='(A,1PE10.3)')'QTABS(JORI,JR  ,JW+1,IC)=',
     &                                QTABS(JORI,JR,JW+1,IC)
            WRITE(0,FMT='(A,1PE10.3)')'QTABS(JORI,JR+1,JW+1,IC)=',
     &                                QTABS(JORI,JR+1,JW+1,IC)
         ENDDO
      ENDIF
!*** diagnostic
!      write(0,*)'readqlib_v19 ckpt 23, wave=',wave
!***
      RETURN
 9000 CONTINUE
      WRITE(0,FMT='(A,I4,A,I4,A,/,A)')
     &   'Fatal error for ICASE=',ICASE,' ISHAP=',ISHAP,
     &   ' :Unable to open file=',PATH//FILENAME(ICASE,ISHAP)
      STOP
 9002 CONTINUE
      WRITE(0,FMT='(A)')'readqlib_v19 ckpt 99: fatal error: stop'
      STOP
 9100 FORMAT(3I4)
 9200 FORMAT(1P8E11.4)
      END
