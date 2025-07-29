      SUBROUTINE STARBURST_RADIATION_FIELD(AGE,UMMP,WAVELENGTH,RADFLD)
      IMPLICIT NONE

! arguments

      DOUBLE PRECISION AGE,RADFLD,UMMP,WAVELENGTH

! local variables

      INTEGER NWAV
      PARAMETER(NWAV=1100)

      INTEGER J,JAGE,JD,JDEV,NSKIP
      LOGICAL INIT
      CHARACTER
     &   DIRNAME*37,
     &   FILENAME*87,
     &   FNAME*50

      DOUBLE PRECISION AGE_STO,KPC,PI,X
      DOUBLE PRECISION
     &   AGE_MOD(1:36),
     &   DUM(1:36),
     &   U(1:NWAV),
     &   USB(1:36),
     &   WAVE(1:NWAV)
      SAVE INIT,JAGE
      SAVE AGE_STO,U,USB,WAVE
      DATA INIT/.TRUE./

!-----------------------------------------------------------------------
! subroutine STARBURST_RADIATION_FIELD
! given
!     WAVELENGTH = wavelength (cm)
!     AGE        = starburst age (yr)

! returns
!     RADFLD = c*u_lambda (erg cm-3 s-1)
!              for 1e6 Msol starburst
!              at a distance of 1 kpc
!     uses results form Starburst99 for Z/Zsol=0.2

!    h = heating rate for 1e-5 cm astrosil grain
!        [MMP83 ISRF with U=1 has h=1.328e-12 erg s-1]

! age      u(1kpc)      L        L      L/M        h     u_mmp
! (Myr)    erg cm-3   erg/s    Lsol   Lsol/Msol erg s-1  
!   1     1.192e-12            1.119e9  1119    1.087e-11  8.185
!   2     1.404e-12            1.317e9  1317    1.310e-11  9.864
!   3     1.814e-12  6.513e42  1.702e9  1702    1.739e-11 13.09
!   4     1.223e-12            1.148e9  1148    1.169e-11  8.803
!   5     9.922e-13  3.562e42  9.311e8   931.1  9.164e-12  6.901
!   6     8.067e-13            7.570e8   757.0  7.099e-12  5.346
!   7     6.534e-13            6.131e8   613.1  5.262e-12  3.962
!   8     5.231e-13            4.909e8   490.9  4.526e-12  3.408
!   9     4.476e-13            4.200e8   420.0  3.917e-12  2.950
!  10     3.864e-13  1.387e42  3.626e8   362.6  3.409e-12  2.567
!  11     3.409e-13            3.199e8   319.9  2.969e-12  2.236
!  12     3.028e-13                             2.580e-12  1.943
!  13     2.734e-13                             2.279e-12  1.716
!  14     2.436e-13                             2.003e-12  1.508
!  15     2.193e-13                             1.809e-12  1.362
!  16     2.036e-13                             1.680e-12  1.265
!  17     1.885e-13                             1.562e-12  1.149
!  18     1.737e-13                             1.448e-12  1.090
!  19     1.630e-13                             1.380e-12  1.039
!  20     1.472e-13                             1.229e-12  0.9255
!  30     9.286e-14                             6.629e-13  0.4992
!  40     6.640e-14                             4.467e-13  0.3364
!  50     5.197e-14                             3.367e-13  0.2535
!  60     4.289e-14                             2.361e-13  0.1778
!  70     3.722e-14                             2.213e-13  0.1666
!  80     3.288e-14                             1.871e-13  0.1409
!  90     2.863e-14                             1.533e-13  0.1154
! 100     2.617e-14                             1.377e-13  0.1037
! 200     1.267e-14                             5.095e-14  0.03837
! 300     8.935e-15                             3.012e-14  0.02268
! 400     6.584e-15                             1.874e-14  0.01411
! 500     5.382e-15                             1.290e-14  0.009714
! 600     4.538e-15                             9.389e-15  0.007070
! 700     3.996e-15                             7.327e-15  0.005517
! 800     3.527e-15                             5.839e-15  0.004397
! 900     3.142e-15                             4.741e-15  0.003570

! B.T.Draine, Princeton University
! history
! 2017.06.12 (BTD) first written
! end history
!-----------------------------------------------------------------------
      DIRNAME='/u/draine/work/lib/astro/starburst99/'
      FNAME='Starburst99_0.2Zsun_instantaneous_starsnebular.dat'
      FILENAME=DIRNAME//FNAME

 1000 CONTINUE

      IF(INIT)THEN
         PI=4.*ATAN(1.D0)
         KPC=3.086E21
         JDEV=17
         OPEN(UNIT=JDEV,FILE=FILENAME)
         JAGE=0.
         DO J=1,20           ! 1,2,...,20
            AGE_MOD(J)=J*1E6
            IF(ABS(AGE/AGE_MOD(J)-1.).LE.0.01)JAGE=J
         ENDDO

         DO J=21,28          ! 30,40,50,...,100
            AGE_MOD(J)=(J-18)*10*1E6
            IF(ABS(AGE/AGE_MOD(J)-1.).LE.0.01)JAGE=J
         ENDDO
         DO J=29,36          ! 200,300,...,900
            AGE_MOD(J)=(J-27)*100*1E6
            IF(ABS(AGE/AGE_MOD(J)-1.).LE.0.01)JAGE=J
         ENDDO

! NSKIP = number of wavelengths to skip (wave < 912 A).

         NSKIP=129
         DO J=1,NSKIP
            READ(JDEV,*)
         ENDDO

         JAGE=0
         DO J=1,36
            IF(ABS(AGE/AGE_MOD(J)-1.).LT.0.01)JAGE=J
         ENDDO
         AGE_STO=AGE

         IF(JAGE.EQ.0)THEN
            write(0,fmt='(a,1pe10.3)')
     &         'starburst_radiation_field ckpt 0 age=',age
            WRITE(0,FMT='(A)')
     &         'starburst_radiation_field fatal error: stop'
            STOP
         ENDIF

! jage= 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20 : 1e6-20e6 yr
         USB(1)=8.185E0
         USB(2)=9.864E0
         USB(3)=13.09E0
         USB(4)=8.803E0
         USB(5)=6.901E0
         USB(6)=5.346E0
         USB(7)=3.962E0
         USB(8)=3.408E0
         USB(9)=2.950E0
         USB(10)=2.567E0
         USB(11)=2.236E0
         USB(12)=1.943E0
         USB(13)=1.716E0
         USB(14)=1.508E0
         USB(15)=1.362E0
         USB(16)=1.265E0
         USB(17)=1.149E0
         USB(18)=1.090E0
         USB(19)=1.039E0
         USB(20)=9.255E-1
         USB(21)=4.992E-1
         USB(22)=3.364E-1
         USB(23)=2.525E-1
         USB(24)=1.778E-1
         USB(25)=1.666E-1
         USB(26)=1.409E-1
         USB(27)=1.154E-1
         USB(28)=1.037E-1
         USB(29)=3.837E-2
         USB(30)=2.268E-2
         USB(31)=1.411E-2
         USB(32)=9.714E-3
         USB(33)=7.070E-3
         USB(34)=5.517E-3
         USB(35)=4.397E-3
         USB(36)=3.570E-3

!*** diagnostic
!         write(0,fmt='(a,1pe10.3)')
!     &      'starburst_radiation_field ckpt 1, age_sto=',age_sto
!         write(0,fmt='(a,i6)')' jage=',jage
!***
            
! NWAV = number of wavelengths to use from file
         DO J=1,NWAV
            READ(JDEV,*)WAVE(J),(DUM(JD),JD=1,36)

!*** diagnostic
!            write(0,fmt='(a,i5,a,1pe10.3)')
!     &         'starburst_radiation_field ckpt 2 j=',j,' wave=',wave(j)
!***
            
! DUM = log10[(dL/dlambda)/(erg s-1 A-1)]
! compute d(flux)/dlambda (erg s-1 cm-3)

            U(J)=1.E8*10.**DUM(JAGE)/(4.*PI*KPC**2)
            IF(UMMP.GT.0.)U(J)=U(J)*UMMP/USB(JAGE)
!*** diagnostic
!            write(0,fmt='(a,i4,a,1pe10.3,a,e10.3)')
!     &         'starburst_radiation_field ckpt 3,j=',j,' dum(jage)=',
!     &         dum(jage),' u(j)=',u(j)
!***

! convert wavelength from A to cm

            WAVE(J)=1.E-8*WAVE(J)

         ENDDO
         CLOSE(JDEV)
         INIT=.FALSE.
      ENDIF
      IF(AGE.NE.AGE_STO)THEN
         INIT=.TRUE.
         GOTO 1000
      ENDIF

      IF(WAVELENGTH.LT.9.11D-6)THEN
         RADFLD=0.
         RETURN
      ELSE

         IF(WAVELENGTH.GT.WAVE(NWAV))THEN
            RADFLD=U(NWAV)*(WAVE(NWAV)/WAVELENGTH)**4
!*** diagnostic
!            write(0,fmt='(a,1pe10.3)')
!     &         'starburst_radiation_field ckpt 3, radfld=',radfld
!***
            RETURN
         ENDIF

! interpolate

         DO J=2,NWAV
            IF(WAVELENGTH.LE.WAVE(J))THEN
               X=(WAVELENGTH-WAVE(J-1))/(WAVE(J)-WAVE(J-1))
               RADFLD=(1.-X)*U(J-1)+X*U(J)
               RETURN
!*** diagnostic
!            write(0,fmt='(a,1pe10.3)')
!     &         'starburst_radiation_field ckpt 4, radfld=',radfld
!            write(0,fmt='(a,1pe10.3)')' wavelength=',wavelength
!            write(0,fmt='(a,i5)')' j=',j
!            write(0,fmt='(a,1pe10.3)')' wave(j-1)=',wave(j-1)
!            write(0,fmt='(a,1pe10.3)')' wave(j)=',wave(j)
!            write(0,fmt='(a,1pe10.3)')' u(j-1)=',u(j-1)
!            write(0,fmt='(a,1pe10.3)')' u(j)=',u(j)
!            write(0,fmt='(a,1pe10.3)')' wave(j)=',wave(j)
!            write(0,fmt='(a,1pe10.3)')' x=',x
!***
            ENDIF
         ENDDO  ! enddo j=2,nwav
            
      ENDIF

!*** diagnostic
!      write(0,fmt='(a)')'starburst_radiation_field : fatal error'
!      write(0,fmt='(a,1pe10.3)')'age=',age
!      write(0,fmt='(a,1pe10.3)')'wavelength=',wavelength
!      stop
!***
      END
