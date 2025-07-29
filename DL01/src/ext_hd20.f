      SUBROUTINE EXT_HD20(WAVEUM,TAUNH)
      IMPLICIT NONE

! Arguments:

      REAL TAUNH,WAVEUM

! Local variables:
! Parameters:
      INTEGER NDIB,NTABMX
      PARAMETER(NDIB=81,NTABMX=10000)

      CHARACTER*43 PATH
      CHARACTER*20 FNAME
      CHARACTER*63 FILENAME
      LOGICAL INIT
      INTEGER I,IUNIT,J,NTAB
      REAL DLOGW,X,Y
      REAL DIBPAR(3,NDIB),XT(1:NTABMX),YT(1:NTABMX)

      DATA INIT/.FALSE./

! DIBPAR(1,J)=wavelength (A) for J-th DIB
! DIBPAR(2,J)=FWHM (A)       for J-th DIB
! DIBPAR(3,J)=EW/E(B-V) (A)  for J-th DIB

      DATA ((DIBPAR(I,J),I=1,3),J=1,40)/
     &   4176.47, 23.33, 0.427,
     &   4428.88, 12.33, 2.231,
     &   4501.80,  2.51, 0.195,
     &   4727.06,  3.92, 0.087,
     &   4761.67, 25.33, 0.623,
     &   4762.57,  2.10, 0.079,
     &   4780.09,  1.45, 0.057,
     &   4824.00, 32.47, 0.452,
     &   4881.83, 19.67, 0.611,
     &   4969.67, 33.70, 0.501,
     &   5109.70, 11.78, 0.274,
     &   5362.14,  5.95, 0.082,
     &   5418.40,  8.07, 0.080,
     &   5449.63, 13.40, 0.244,
     &   5456.00, 27.05, 0.134,
     &   5487.43,  5.36, 0.121,
     &   5487.49,  1.79, 0.069,
     &   5508.35,  3.54, 0.132,
     &   5524.89,  3.73, 0.054,
     &   5535.68, 12.92, 0.073,
     &   5537.00, 22.57, 0.285,
     &   5704.75,  6.64, 0.081,
     &   5705.13,  2.15, 0.096,
     &   5779.48, 15.50, 0.647,
     &   5780.59,  2.07, 0.579,
     &   5795.23,  4.09, 0.117,
     &   5797.11,  0.97, 0.132,
     &   5844.19,  3.04, 0.077,
     &   6010.58,  3.46, 0.141,
     &   6045.27, 14.23, 0.189,
     &   6177.27, 23.00, 0.773,
     &   6196.19,  0.65, 0.061,
     &   6203.19,  1.21, 0.107,
     &   6204.33,  3.93, 0.189,
     &   6207.83, 11.70, 0.136,
     &   6270.06,  1.03, 0.076,
     &   6270.36,  3.56, 0.061,
     &   6281.07,  8.47, 1.237,
     &   6284.31,  2.58, 0.618,
     &   6308.93,  2.37, 0.058/
      DATA ((DIBPAR(I,J),I=1,3),J=41,81)/
     &   6315.00, 21.23, 0.352,
     &   6317.58,  2.06, 0.058,
     &   6359.53, 37.33, 0.536,
     &   6379.27,  0.79, 0.078,
     &   6413.50,  8.07, 0.085,
     &   6451.60, 25.40, 0.403,
     &   6494.17, 11.17, 0.200,
     &   6532.10, 17.20, 0.664,
     &   6591.40,  5.55, 0.087,
     &   6613.72,  1.14, 0.231,
     &   6660.64,  0.84, 0.051,
     &   6919.25,  0.96, 0.053,
     &   6939.00, 21.27, 0.396,
     &   6993.18,  0.96, 0.116,
     &   7223.13,  5.36, 0.083,
     &   7224.18,  1.07, 0.259,
     &   7274.50,  5.06, 0.058,
     &   7334.33,  1.24, 0.060,
     &   7357.20, 28.23, 0.227,
     &   7432.07, 22.33, 0.549,
     &   7562.24,  1.78, 0.087,
     &   7569.70,  5.50, 0.216,
     &   7651.37,  2.72, 0.083,
     &   7686.46,  3.17, 0.073,
     &   7705.90,  3.53, 0.081,
     &   7709.67, 33.53, 0.444,
     &   7721.03,  1.81, 0.067,
     &   7748.18,  4.53, 0.072,
     &   7782.25,  3.64, 0.072,
     &   7843.87,  4.80, 0.083,
     &   7927.80, 15.00, 0.428,
     &   7987.89,  2.90, 0.051,
     &   8038.48,  3.20, 0.084,
     &   8350.79,  1.73, 0.073,
     &   8621.11,  1.86, 0.125,
     &   8621.23,  5.59, 0.272,
     &   8648.28,  4.17, 0.241,
     &   9577.0,   4.1,  0.398,
     &   9632.0,   4.0,  0.573,
     &   11797.0,  2.7,  0.13,
     &   13175.0,  4.0,  0.36/

      SAVE DIBPAR,DLOGW,NTAB,XT,YT

!***********************************************************************
! Given:
!       WAVEUM    = wavelength (um) 
! Returns:
!       TAUNH = tau_ext/N_H (cm^2/H)
!
! Optical DIBs:
!
! We include 77 DIB features between 8680-3800 A with EW/E(B-V)>0.05A
! classified as "certain" (69) or "probable" (8) by Jenniksens & 
! Desert 1994, with central wavelength, width, and EW/E(B-V)
! taken from the catalogue of Jenniksens & Desert 1994
! plus bands at 9577A, 9632A, 1.18um, 1.32um
!
! Diffuse interstellar bands are added assuming Drude profiles
!
! 9577A, 9632 band params from updated catalog of Jenniskens
!       http://www-space.arc.nasa.gov/~leonid/DIBcatalog.html  [2003.04.01]
!       based on Foing & Ehrenfreund 1994, Nature 369, 296
! 
! 1.1797um band from Joblin etal 1990 (Nature 346, 729)
! 1.3175um band from Joblin etal 1990 (Nature 346, 729).
!
! Current version of code asssumes that tau(DIB)/N_H is independent of 
! environment
! References:
!
! Hensley, B.S., \& Draine, B.T. 2020, ApJ, submitted
!
! B.T. Draine, Princeton Univ. Obs., 2020.09.10
! History:
! 20.09.10 (BTD) created
! end history
!***********************************************************************

!*** diagnostic
!      write(0,*)'ext_HD20 ckpt 0, wave=',wave
!***
      IF(.NOT.INIT)THEN
! Initialization:
!         PATH= '/u/draine/work/papers/hensley+draine_2020b/'
         PATH= '/u/draine/work/papers/hensley+draine_2021a/'
         FNAME='tau_nh_31aug2020.dat'
         FILENAME=PATH//FNAME
         IUNIT=37
         OPEN(UNIT=IUNIT,FILE=FILENAME,STATUS='old')
         READ(IUNIT,*)
         J=0
 1000    READ(IUNIT,*,END=2000,ERR=2000)X,Y
         J=J+1
         XT(J)=X
         YT(J)=Y
         NTAB=J
         GOTO 1000
 2000    CLOSE(IUNIT)
         INIT=.TRUE.
         DLOGW=LOG(XT(NTAB)/XT(1))/(NTAB-1)
      ENDIF

      
! ************ Initialization is complete ***************

! Sanity check.  We will allow use of ext_HD20.f only for E < 13.6eV

      IF(WAVEUM.LT.0.0911)THEN
         WRITE(0,*)'ERROR: ext_HD20.f not ready for lambda < 911 A!'
         STOP
      ELSE

! we assume that XT are logarithmically spaced

         X=LOG(WAVEUM/XT(1))/LOG(XT(NTAB)/XT(1))
         J=INT(1.+(NTAB-1)*X)
         IF(J.LE.0)J=1
         IF(J.GE.NTAB)J=NTAB-1
         X=LOG(WAVEUM/XT(J))/DLOGW
!*** sanity check
         IF(X.LT.0.OR.X.GT.1)THEN
            WRITE(0,FMT='(A,1PE10.3,A,1PE10.3)')
     &         'ext_hd20 ckpt 9 fatal error: x=',X,' WAVEUM=',WAVEUM
            STOP
         ENDIF
!***
         TAUNH=EXP((1.-X)*LOG(YT(J))+X*LOG(YT(J+1)))
      ENDIF

      RETURN
      END
