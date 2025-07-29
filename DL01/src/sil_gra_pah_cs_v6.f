      SUBROUTINE SIL_GRA_PAH_CRS_SCT(ICASE,BOVERA,ACM,DFRAC,
     &                               WAVELENGTHCM,CABS,CSCA)
      IMPLICIT NONE
      CHARACTER COMPOSITION*3
      INTEGER ICASE,ISHAP
      DOUBLE PRECISION ACM,BOVERA,CABS,CSCA,DFRAC,WAVELENGTHCM

! local variables:

      INTEGER ICGREMT2
      REAL AUM,BOVERA_SP,G,G2,QABS,QEXT,QSCA,TGR,WAVEUM
      REAL 
     &   QABSV(3),
     &   QEXTV(3),
     &   QSCAV(3)

! ------------- sil_gra_pah_cs_v6.f ----------------------------
! input:    
!    ICASE = designator of grain material (to be passed to READQLIB)
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
!              21 for D16 graphite-PAH0, c || symm axis a, a_t=50A
!              22 for D16 graphite-PAH+, EMT, a_t=20A
!              23 for D16 graphite-PAH0, EMT, a_t=20A
!              24 for D16 graphite-PAH+, EMT, a_t=50A
!              25 for D16 graphite-PAH0, EMT, a_t=50A
!              26 for D16 graphite, random c axis
!              27 for D16 graphite, c || a
!              28 for D16 graphite, EMT (MG1: MG with matrix=E||c)
!          29-118 for DH19 Astrodust
!                 ICASE=29+IDX
!                 IDX= JPO*6 + JFE
!                 JPO=0-9 for porosity = 0, 0.1, 0.2, ..., 0.9  [10 cases]
!                 JFE=0-5 for f_Fe = 0, 0.1, 0.2, ... 0.5   [ 6 cases]
!             119 for D16 graphite, EMT (MG2: MG with matrix=Eperpc)

!       A           =  equal-volume radius (cm)
!       DFRAC        = D/(H+D) ratio in PAHs;
!       WAVELENGTHCM = wavelength (cm)

! output:
!       CABS = absorption cross section (cm^2)
!       CSCA = scattering cross section (cm^2)
!
! 
! ***************************************************************
! Calculate absorption cross section of graphite/silicate.    -
! Created by Aigen Li, Princeton Univ. Obs., 1999.10.11
! History:
! 1999.10.11 (AL)  First written
! 1999.11.29 (AL)  include both PAH and PAH+
! 2000.03.08 (AL)  add deuterated PAHs
! 2005.10.29 (BTD) modify for compatibility with new version
!                  of QCOMP
!                  cosmetic modification
! 2006.06.05 (BTD) modify for consistency with new deuteration
!                  options in QCOMP [D/(H+D)=0.25 or 0.5]
! 2006.06.13 (BTD) changed DTOH to DFRAC to avoid confusion
! 2007.12.10 (BTD) modified to recognize amc option
! 2007.12.23 (BTD) changed
!                  ICOMP=46->40
!                        47->41
!                  for consisistency with new qcomp_v4.f
! 2007.12.24 (BTD) modifiy to use new 2-component model for
!                  graphite E||c for both "gra" and "PAH"
!                  options
! 2011.02.19 (BTD) v3
!                  * modified to support composition 'fem'
! 2011.02.23 (BTD) * further modifications
! 2014.02.17 (BTD) v4
!                  * now simply an interface to subroutine
!                    READQLIB
! 2017.03.05 (BTD) * add support for ICASE=34,35,36,37,38
! 2017.10.10 (BTD) * replaced ISHAP with BOVERA in arg list 
!                  * replaced ISHAP with BOVERA_SP in arg list of
!                    READQLIB
! 2017.12.25 (BTD) * updated comments to describe ICASE=39-42
! 2019.01.21 (BTD) v5
!                  * add support for ICASE=569
! 2019.08.01 (BTD) v6
!                  * change comments (but not code) for use with
!                    DH19 Astrodust instead of DH18 Astrosil
!                    [need to use readqlib_v18.f]
!     end history
!
! ---------------------------------------------------------------
!*** diagnostic
!      write(0,fmt='(a,i3)')'sil_gra_pah_cs_v6 ckpt 1: icase=',
!     &                icase
!      write(0,fmt='(a,1pe10.3)')' bovera=',bovera
!      write(0,fmt='(a)')' about to call readqlib'
!***
! ---------------------------------------------------------------

! BE CAREFUL: wavelength, aum are single precision!!

! wavelengthcm: cm; wavelength: micron;

      WAVEUM=SNGL(WAVELENGTHCM*1.E4)

      AUM=SNGL(ACM*1.E4)
      BOVERA_SP=SNGL(BOVERA)

! a: size in cm; a_um in micron; 		

      CALL READQLIB(ICASE,BOVERA_SP,AUM,WAVEUM,QABSV,QEXTV,QSCAV)
      QABS = (QABSV(1)+QABSV(2)+QABSV(3))/3.0
      QSCA = (QSCAV(1)+QSCAV(2)+QSCAV(3))/3.0
      CABS=3.141593D0*ACM*ACM*QABS
      CSCA=3.141593D0*ACM*ACM*QSCA

      RETURN
      END
