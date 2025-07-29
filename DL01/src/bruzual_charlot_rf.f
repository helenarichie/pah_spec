      SUBROUTINE BRUZUAL_CHARLOT_RF(AGE,ZBC,UMMP,WAVELENGTH,RADFLD)
      
      IMPLICIT NONE

! arguments

      DOUBLE PRECISION AGE,RADFLD,UMMP,WAVELENGTH,ZBC

! local variables
! NSP = number of different stellar population models

      INTEGER NSP,NWAVA
      PARAMETER(NSP=6,NWAVA=1200)

      INTEGER J,JD,JDEV,JSP,NSKIP,NWAV
      LOGICAL INIT
      CHARACTER
     &   DIRNAME*43,
     &   FILENAME*87
      CHARACTER
     &   FNAME(1:NSP)*50

      DOUBLE PRECISION AGE_STO,FAC,KPC,PI,USUM,X,ZBC_STO
      DOUBLE PRECISION
     &   U(1:NWAVA),
     &   USP(1:NSP),
     &   WAVE(1:NWAVA)
      SAVE INIT,JSP,NWAV
      SAVE AGE_STO,U,USP,WAVE,ZBC_STO
      DATA INIT/.TRUE./

!-----------------------------------------------------------------------
! subroutine BRUZUAL_CHARLOT_RADIATION_FIELD
! given
!     AGE        = starburst age (yr)
!     ZBC        = heavy element mass fraction
!     UMMP       = (dust heating rate) / (heating rate for mMMP RF)
!     WAVELENGTH = wavelength (cm)

! returns
!     RADFLD = c*u_lambda (erg cm-3 s-1)

!     uses Bruzual-Charlot models for single-age stellar population

!    h = heating rate for 1e-5 cm astrosil grain
!        [mMMP ISRF with U=1 has h=1.974e-12 erg s-1 for DH20Ad_P0.20_0.00


! B.T.Draine, Princeton University
! history
! 2019.08.02 (BTD) first written
!                  adapted from starburst_radiation_field.f
! 2020.05.04 (BTD) added new case age=2.86e8 yr 
!                  (from Leslie Hunt)
!                  now have 6 cases
! 2020.05.21 (BTD) modified to begin reading tabulated spectrum
!                  at wavelength just shortward of Lyman limit
! 2020.09.20 (BTD) further adjust intensity to match h for mMMP
!     end history
!-----------------------------------------------------------------------
      DIRNAME='/u/draine/work/lib/astro/bruzual-charlot03/'
!***  diagnostic
!      write(0,fmt='(a)')'bruzual_charlot ckpt 0'
!      write(0,fmt='(a,a)')'  dirname=',dirname
!***
      FNAME(1)='bc03_z0.02_age3.0200e+06.dat'
      FNAME(2)='bc03_z0.02_age1.0000e+07.dat'
      FNAME(3)='bc03_z0.02_age1.0152e+08.dat'
      FNAME(4)='bc03_z0.02_age2.8612e+08.dat'
      FNAME(5)='bc03_z0.02_age1.0152e+09.dat'
      FNAME(6)='bc03_z0.0004_age1.0000e+07.dat'

!***  diagnostic
!      write(0,fmt='(a)')'bruzual_charlot ckpt 1'
!      write(0,fmt='(a,a)')'  fname(5)=',fname(5)
!***
      IF(.NOT.INIT)THEN
         IF(ZBC.NE.ZBC_STO.OR.AGE.NE.AGE_STO)THEN
            INIT=.TRUE.
         ENDIF
      ENDIF

      IF(INIT)THEN
         IF(ABS(ZBC-0.02).LE.1.E-3)THEN
            IF(AGE.GT.1.E6.AND.AGE.LE.5.E6)THEN
               JSP=1   ! Z=0.02, t=3e6
            ELSEIF(AGE.GT.5.E6.AND.AGE.LE.3.E7)THEN
               JSP=2   ! Z=0.02, t=1e7
            ELSEIF(AGE.GT.3.E7.AND.AGE.LE.1.7E8)THEN
               JSP=3   ! Z=0.02, t=1e8
            ELSEIF(AGE.GT.1.7E8.AND.AGE.LE.5.5E8)THEN
               JSP=4   ! Z=0.02, t=3e8
            ELSEIF(AGE.GT.5.5E8.AND.AGE.LE.3.E9)THEN
               JSP=5   ! Z=0.02, t=1e9
            ENDIF
         ELSEIF(ABS(ZBC-0.0004).LE.0.0001)THEN
            JSP=6      ! Z=0.0004, t=1e7
         ELSE
            WRITE(0,FMT='(A,1PE10.3,A,1PE10.3,A)')
     &         'bruzual_charlot_rf ckpt 3: problem: for age=',AGE,
     &         ' and ZBC=',ZBC,' : halt execution'
            STOP
         ENDIF
         
!***
         WRITE(0,FMT='(A,I3)')'bruzual_charlot_rf ckpt 2: jsp=',jsp
!***

         FILENAME=DIRNAME//FNAME(JSP)
         
         AGE_STO=AGE
         ZBC_STO=ZBC
         
         PI=4.*ATAN(1.D0)
         KPC=3.086E21
         JDEV=17
!***
         write(0,fmt='(a,a)')'bruzual_charlot_rf ckpt 4: filename=',
     &      filename
!***
         OPEN(UNIT=JDEV,FILE=FILENAME)

! NSKIP = number of wavelengths to skip (wave < 911 A).

         NSKIP=123
         DO J=1,NSKIP
            READ(JDEV,*)
         ENDDO

! NWAV = number of wavelengths to use from file
! convert wavelength from A to cm

         DO J=1,NWAVA
            READ(JDEV,*,END=1000)WAVE(J),U(J)
            WAVE(J)=1.E-8*WAVE(J)
            NWAV=J
         ENDDO
 1000    CONTINUE

!*** diagnostic
         WRITE(0,fmt='(a,i6)')'bruzual_charlot_rf ckpt 4.5, nwav=',nwav
!***

! USP(JSP) = <C_abs>_JSP / <C_abs>_MMP83
         
!         USP(1)=7.430E0   ! Z=0.02, t=3e6 yr
!         USP(2)=6.032E0   ! Z=0.02, t=1e7 yr
!         USP(3)=3.830E0   ! Z=0.02, t=1e8 yr
!         USP(4)=2.213E0   ! Z=0.02, t=3e8 yr
!         USP(5)=1.039E0   ! Z=0.02, t=1e9 yr
!         USP(6)=7.432E0   ! Z=0.0004, t=1e7 yr

! corrected 2020.09.21 (BTD):

!         USP(1)=5.606E0   ! Z=0.02, t=3e6 yr
!         USP(2)=4.567E0   ! Z=0.02, t=1e7 yr
!         USP(3)=2.901E0   ! Z=0.02, t=1e8 yr
!         USP(4)=1.676E0   ! Z=0.02, t=3e8 yr
!         USP(5)=0.787E0   ! Z=0.02, t=1e9 yr
!         USP(6)=5.623E0   ! Z=0.0004, t=1e7 yr
         
! corrected 2021.02.01 (BTD) for b/a=1.6:

         USP(1)=5.549E0   ! Z=0.02, t=3e6 yr
         USP(2)=4.530E0   ! Z=0.02, t=1e7 yr
         USP(3)=2.900E0   ! Z=0.02, t=1e8 yr
         USP(4)=1.688E0   ! Z=0.02, t=3e8 yr
         USP(5)=0.7944E0  ! Z=0.02, t=1e9 yr
         USP(6)=5.566E0   ! Z=0.0004, t=1e7 yr
         
         USUM=0.5*(U(1)*(WAVE(2)-WAVE(1))+
     &             U(NWAV)*(WAVE(NWAV)-WAVE(NWAV-1)))
         DO J=2,NWAV-1
            USUM=USUM+U(J)*0.5*(WAVE(J+1)-WAVE(J-1))
         ENDDO
         USUM=USUM/2.9979E10
         
! a=1e-5cm P=0.20 b/a=1.6 astrodust grain in mMMP radiation field
! is heated at rate h=1.958e-12 erg s-1
! (a=1e-5cm P=0.20 b/a=0.5 astrodust grain in mMMP radiation field
! is heated at rate h=1.974e-12 erg s-1)
! mMMP radiation field has energy density 1.043e-12 erg cm-3
! compute U(J) = dflux/dlambda (erg s-1 cm-3)
! adjusted to give heating rate for a=1e-5cm b/a=1.6 P=0.20 astrodust grain
! equal to UMMP times the heating rate for the same grain in the
! mMMP radiation field

!         FAC=UMMP*1.276E-12/(USP(JSP)*USUM)
         FAC=UMMP*1.043E-12/(USP(JSP)*USUM)
         DO J=1,NWAV
            U(J)=FAC*U(J)
         ENDDO

         CLOSE(JDEV)
         INIT=.FALSE.
      ENDIF   ! endif(INIT)

      IF(WAVELENGTH.LT.9.11D-6)THEN
         RADFLD=0.
         RETURN
      ELSE

         IF(WAVELENGTH.GT.WAVE(NWAV))THEN
            RADFLD=U(NWAV)*(WAVE(NWAV)/WAVELENGTH)**4
            RETURN
         ENDIF

! interpolate

         DO J=2,NWAV
            IF(WAVELENGTH.LE.WAVE(J))THEN
               X=(WAVELENGTH-WAVE(J-1))/(WAVE(J)-WAVE(J-1))
               RADFLD=(1.-X)*U(J-1)+X*U(J)
!*** diagnostic
               if(radfld.lt.0.)then
                  write(0,fmt='(a,1pe10.3)')
     &               'bruzual-charlot_rf ckpt 1: wavelength=',wavelength
                  write(0,fmt='(a,1pe10.3)')'radfld=',radfld
                  write(0,fmt='(a,i4,a,1pe10.3,a,e10.3)')
     &               'j-1=',j-1,' wave(j-1)=',wave(j-1),
     &               ' u(j-1)=',u(j-1)
                  write(0,fmt='(a,i4,a,1pe10.3,a,e10.3)')
     &               'j  =',j,  ' wave(j) = ',wave(j),' u(j) =',u(j)
                  write(0,fmt='(a)')'fatal error'
                  stop
               endif
!***
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
      RETURN
      END
