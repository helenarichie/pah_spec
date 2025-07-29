      SUBROUTINE KURUCZ_RF(TEFF,UMMP,WAVELENGTH,RADFLD)
      
      IMPLICIT NONE

! arguments

      DOUBLE PRECISION RADFLD,TEFF,UMMP,WAVELENGTH

! local variables
! NSP = number of different stellar population models

      INTEGER NSP,NWAVA
      PARAMETER(NSP=11,NWAVA=1000)

      INTEGER J,JD,JDEV,JSP,NSKIP,NWAV
      LOGICAL INIT
      CHARACTER
     &   DIRNAME*32,
     &   FILENAME*87
      CHARACTER
     &   FNAME(1:NSP)*50

      DOUBLE PRECISION FAC,KPC,PI,TEFF_STO,USUM,X
      DOUBLE PRECISION
     &   U(1:NWAVA),
     &   USP(1:NSP),
     &   WAVE(1:NWAVA)
      SAVE INIT,JSP,NWAV
      SAVE TEFF_STO,U,USP,WAVE
      DATA INIT/.TRUE./

!-----------------------------------------------------------------------
! subroutine BRUZUAL_CHARLOT_RADIATION_FIELD
! given
!     TEFF       = effective temperature
!     UMMP       = (dust heating rate) / (heating rate for MMP83 RF)
!     WAVELENGTH = wavelength (cm)

! returns
!     RADFLD = c*u_lambda (erg cm-3 s-1)

!     uses Kurucz stellar atmospheres
!          for solar metallicity and main-sequence gravity

!    h = heating rate for 1e-5 cm astrosil grain
!        [MMP83 ISRF with U=1 has h=1.328e-12 erg s-1 for DHsil...
!                                   1.625e-12 erg s-1 for DH19Ad_P0.20_0.00
!                                   1.829e-12 erg s-1 for DH20Ad_P0.20_0.00


! B.T.Draine, Princeton University
! history
! 2020.08.11 (BTD) adapted from bruzual_charlot_rf.f
!                  adapted from starburst_radiation_field.f
!     end history
!-----------------------------------------------------------------------
      DIRNAME='/u/draine/work/lib/astro/kurucz/'   ! 32 characters
!***  diagnostic
!      write(0,fmt='(a)')'kurucz_rf ckpt 0'
!      write(0,fmt='(a,a)')'  dirname=',dirname
!***
      FNAME(1)='4piJlambda_T3.0E3_g4.0_z0.0.dat'
      FNAME(2)='4piJlambda_T3.5E3_g4.0_z0.0.dat'
      FNAME(3)='4piJlambda_T4.5E3_g4.0_z0.0.dat'
      FNAME(4)='4piJlambda_T5.0E3_g4.0_z0.0.dat'
      FNAME(5)='4piJlambda_T5780_g4.4_z0.0.dat'
      FNAME(6)='4piJlambda_T6.0E3_g4.0_z0.0.dat'
      FNAME(7)='4piJlambda_T8.0E3_g4.0_z0.0.dat'
      FNAME(8)='4piJlambda_T1.0E4_g4.0_z0.0.dat'
      FNAME(9)='4piJlambda_T1.5E4_g4.0_z0.0.dat'
      FNAME(10)='4piJlambda_T2.2E4_g4.0_z0.0.dat'
      FNAME(11)='4piJlambda_T3.0E4_g4.0_z0.0.dat'

!***  diagnostic
!      write(0,fmt='(a)')'kurucz_rf ckpt 1'
!      write(0,fmt='(a,a)')'  fname(5)=',fname(5)
!***
      IF(.NOT.INIT)THEN
         IF(TEFF.NE.TEFF_STO)THEN
            INIT=.TRUE.
         ENDIF
      ENDIF

      IF(INIT)THEN   ! need to initialize
         IF(TEFF.GT.2.90E3.AND.TEFF.LE.3.25E3)THEN
            JSP=1   ! Teff=3000
         ELSEIF(TEFF.GT.3.25E3.AND.TEFF.LE.4.00E3)THEN
            JSP=2   ! Teff=3500
         ELSEIF(TEFF.GT.4.00E3.AND.TEFF.LE.4.75E3)THEN
            JSP=3   ! Teff=4500
         ELSEIF(TEFF.GT.4.75E3.AND.TEFF.LE.5.50E3)THEN
            JSP=4   ! Teff=5000
         ELSEIF(TEFF.GT.5.50E3.AND.TEFF.LE.5.85E3)THEN
            JSP=5   ! Teff=5780
         ELSEIF(TEFF.GT.5.85E3.AND.TEFF.LE.7.00E3)THEN
            JSP=6   ! Teff=6000
         ELSEIF(TEFF.GT.7.00E3.AND.TEFF.LE.9.00E3)THEN
            JSP=7   ! Teff=8000
         ELSEIF(TEFF.GT.9.00E3.AND.TEFF.LE.12.5E3)THEN
            JSP=8   ! Teff=10000
         ELSEIF(TEFF.GT.12.5E3.AND.TEFF.LE.18.5E3)THEN
            JSP=9   ! Teff=15000
         ELSEIF(TEFF.GT.18.5E3.AND.TEFF.LE.26.0E3)THEN
            JSP=10   ! Teff=22000
         ELSEIF(TEFF.GT.26.0E3.AND.TEFF.LE.32.0E3)THEN
            JSP=11   ! Teff=30000
         ELSE
            WRITE(0,FMT='(A,1PE10.3,A)')
     &         'kurucz_rf ckpt 3: problem: for Teff=',TEFF,
     &         ' : halt execution'
            STOP
         ENDIF
 
!***
         WRITE(0,FMT='(A,I3)')'kurucz_rf ckpt 2: jsp=',jsp
!***

         FILENAME=DIRNAME//FNAME(JSP)
         
         TEFF_STO=TEFF
         
         PI=4.*ATAN(1.D0)

         JDEV=17
!***
         write(0,fmt='(a,a)')'kurucz_rf ckpt 4: filename=',
     &      filename
!***
         OPEN(UNIT=JDEV,FILE=FILENAME)

! NSKIP = number of wavelengths to skip (wave < 911 A).

         NSKIP=6
         DO J=1,NSKIP
            READ(JDEV,*)
         ENDDO

! NWAV = number of wavelengths to use from file
! convert wavelength from um to cm
! input: U(J)=4*pi*J_lambda  (erg cm-3 s-1)
         DO J=1,NWAVA
            READ(JDEV,*,END=1000)WAVE(J),U(J)
            WAVE(J)=1.E-4*WAVE(J)
            NWAV=J
         ENDDO
 1000    CONTINUE

!*** diagnostic
         WRITE(0,fmt='(a,i6)')'kurucz_rf ckpt 5, nwav=',nwav
!***

! Precalculated:
! USP(JSP) = <C_abs>_JSP / <C_abs>_MMP83
         
         USP(1)=0.519     ! Teff=3000 ok
         USP(2)=0.5517    ! Teff=3500 ok
         USP(3)=0.7496    ! Teff=4500 ok
         USP(4)=0.870     ! Teff=5000 ok
         USP(5)=1.125     ! Teff=5780 ok
         USP(6)=1.131     ! Teff=6000
         USP(7)=1.722     ! Teff=8000
         USP(8)=2.974     ! Teff=10000
         USP(9)=5.385     ! Teff=15000
         USP(10)=6.840    ! Teff=22000
         USP(11)=7.408    ! Teff=30000
         
         USUM=0.5*(U(1)*(WAVE(2)-WAVE(1))+
     &             U(NWAV)*(WAVE(NWAV)-WAVE(NWAV-1)))
         DO J=2,NWAV-1
            USUM=USUM+U(J)*0.5*(WAVE(J+1)-WAVE(J-1))
         ENDDO
         USUM=USUM/2.9979E10
         
! a=1e-5cm P=0.20 b/a=0.5 astrodust grain in MMP83 radiation field
! is heated at rate h=1.829e-12 erg s-1
! MMP83 radiation field has energy density 1.276e-12 erg cm-3
! compute U(J) = dflux/dlambda (erg s-1 cm-3)
! adjusted to give heating rate for a=1e-5cm astrodust grain
! equal to UMMP times the heating rate for the same grain in the
! MMP83 radiation field

         FAC=UMMP*1.276E-12/(USP(JSP)*USUM)
         DO J=1,NWAV
            U(J)=FAC*U(J)
         ENDDO

         CLOSE(JDEV)
!*** diagnostic
!         WRITE(0,fmt='(a)')'kurucz_rf ckpt 6'
!***
         INIT=.FALSE.
      ENDIF   ! endif(INIT)

      IF(WAVELENGTH.LT.9.11D-6)THEN
         RADFLD=0.
         RETURN
      ELSE

! Rayleigh-Jeans limit: u_lambda propto lambda^{-4}

         IF(WAVELENGTH.GT.WAVE(NWAV))THEN
            RADFLD=U(NWAV)*(WAVE(NWAV)/WAVELENGTH)**4
         ELSE

! interpolate

            DO J=2,NWAV

! tabulated spectrum is in order of increasing wavelength
! scan through list until entries bracket wavelength

!*** diagnostic
!         WRITE(0,fmt='(a,i6)')'kurucz_rf ckpt 7, j=',j
!         write(0,fmt='(a,1pe10.3,a,1pe10.3)')
!     &      '   wavelength=',wavelength,' wave(j)=',wave(j)
!***
               IF(WAVELENGTH.LE.WAVE(J))THEN
!*** diagnostic
!                  write(0,fmt='(a)')
!     &               'kurucz_rf ckpt 8'
!                  write(0,fmt='(a,1pe10.3)')
!     &               '   wavelength=',wavelength
!                  write(0,fmt='(a,1pe10.3)')
!     &               '   wave(j)=',wave(j)
!                  write(0,fmt='(a,1pe10.3)')
!     &               '   wave(j-1)=',wave(j-1)
!***
                  X=(WAVELENGTH-WAVE(J-1))/(WAVE(J)-WAVE(J-1))
!*** diagnostic
!                  write(0,fmt='(a,1pe10.3,a,1pe10.3,a,1pe10.3)')
!     &               'kurucz_rf ckpt 8.1, x=',x,' u(j-1)=',u(j-1),
!     &               ' u(j)=',u(j)
!***
                  RADFLD=(1.-X)*U(J-1)+X*U(J)
!*** diagnostic
!                  write(0,fmt='(a,1pe10.3,a,1pe10.3)')
!     &               'kurucz_rf ckpt 9, x=',x,' radfld=',radfld
!*** diagnostic
                  if(radfld.lt.0.)then
                     write(0,fmt='(a,1pe10.3)')
     &                  'kurucz_rf ckpt 10: wavelength=',wavelength
                     write(0,fmt='(a,1pe10.3)')'radfld=',radfld
                     write(0,fmt='(a,i4,a,1pe10.3,a,e10.3)')
     &                  'j-1=',j-1,' wave(j-1)=',wave(j-1),
     &                  ' u(j-1)=',u(j-1)
                     write(0,fmt='(a,i4,a,1pe10.3,a,e10.3)')
     &                  'j  =',j,  ' wave(j) = ',wave(j),' u(j) =',u(j)
                     write(0,fmt='(a)')'fatal error'
                     stop
                  endif
!***    
                  RETURN
               ENDIF
            ENDDO  ! enddo j=2,nwav

!*** diagnostic: should not get here
!            write(0,fmt='(a,a)')
!     &         'kurucz_rf ckpt 11, PROBLEM: Teff=',TEFF
!            write(0,fmt='(a,1pe10.3)')' wavelength=',wavelength
            STOP
!***
         ENDIF
            
      ENDIF

!*** diagnostic
      write(0,fmt='(a)')'kurucz_rf : fatal error'
      write(0,fmt='(a,1pe10.3)')'Teff=',TEFF
      write(0,fmt='(a,1pe10.3)')'wavelength=',wavelength
      stop
!***
      RETURN
      END
