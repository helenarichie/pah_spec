      SUBROUTINE SILPROF_GC(WAVE,SILPROF)
      IMPLICIT NONE

! --------------------- silprof_gc_v2 ---------------------------------
! arguments

      REAL SILPROF,WAVE

! local variables

      LOGICAL INIT
      INTEGER J,JHI,JLO,JT
      REAL GAMMA2,X
      REAL 
     &   PROFTAB(27),
     &   WAVTAB(27)

      SAVE INIT,PROFTAB,WAVTAB

!=======================================================================
! Subroutine SILPROF_GC
! Given:
!       WAVE = wavelength (um)
! Returns:
!       SILPROF = (Delta tau)_lambda/(Delta tau)_9.7
!                 contributed by silicate profile
!                 using profile observed by 
!                 Kemper, Vriend, and Tielens (2004)
!                 toward Galactic Center (ApJ 609, 826)
! history:
! 05.03.25 (BTD) first written for use by extf99.f
! 10.01.17 (BTD) v2
!                * add contribution of 18um drude profile for
!                  lambda < 12.7um for continuity
!                * modified file GC_sil_prof.dat to allow
!                  for contribution of 18um profile over 8-12.6um
!                  range
! 14.07.04 (BTD) v2
!                * peel off from ext_f99 as separate routine
!                  in order to use with other code
! end history
!=======================================================================
      DATA INIT/.TRUE./

!*** diagnostic
!      write(0,*)'silprof_gc_v2 ckpt 0: wave=',wave
!***
      IF(INIT)THEN
!*** diagnostic
!         write(0,*)'silprof_gc_v2 ckpt 1'
!***
         OPEN(UNIT=33,FILE='/u/draine/work/opt/dat/GC_sil_prof.dat')
         DO J=1,3
            READ(33,*)
         ENDDO
         DO J=1,27
            READ(33,*)WAVTAB(J),PROFTAB(J)
         ENDDO
         CLOSE(33)
         INIT=.FALSE.
      ENDIF

! On short wavelength side, assume that profile declines exponentially
! with wavelength shortward of 8.0um (factor of 1.5 each 0.2um)

      IF(WAVE.LE.WAVTAB(1))THEN
!*** diagnostic
!         write(0,*)'ext_f99_v2 silprov_gc ckpt 2 (short wave)'
!***
         SILPROF=(PROFTAB(1)+0.0142)*1.5**(5.*(WAVE-WAVTAB(1)))
         RETURN

      ELSEIF(WAVE.GE.WAVTAB(27))THEN

! On long wavelength side, assume drude profile
!     0.4*gamma^2/[(w/18-18/w)**2+gamma^2]
! gamma^2=0.12 gives a reasonable join to the galactic center profile
! FWHM=sqrt(gamma2)*18=6.235um

!*** diagnostic
!         write(0,*)'silprof_gc_v2 ckpt 3 (long wave)'
!***

         GAMMA2=0.12
         SILPROF=0.4*GAMMA2/((WAVE/18.-18./WAVE)**2+GAMMA2)
         RETURN

      ELSE
!*** diagnostic
!         write(0,*)'silprof_gc_v2 ckpt 4'
!***
         JLO=1
         JHI=27
         JT=14
      ENDIF

 1000 CONTINUE

!*** diagnostic
!      write(0,*)'silprof_gc_v2 ckpt 5'
!***

      IF(WAVE.GE.WAVTAB(JT))THEN
         JLO=JT
      ELSE
         JHI=JT
      ENDIF
      IF((JHI-JLO).GT.1)THEN
         JT=(JLO+JHI)/2
         GOTO 1000
      ENDIF

!*** diagnostic
!      write(0,*)'silprof_gc_v2 ckpt 6: jlo,jhi=',jlo,jhi
!***

! have bracketed WAVE.
! do linear interpolation

      X=(WAVE-WAVTAB(JLO))/(WAVTAB(JHI)-WAVTAB(JLO))
      SILPROF=(1.-X)*PROFTAB(JLO)+X*PROFTAB(JHI)

! add contribution from 18um drude profile:

      GAMMA2=0.12
      SILPROF=SILPROF+0.4*GAMMA2/((WAVE/18-18./WAVE)**2+GAMMA2)

!*** diagnostic
!      write(0,*)'wave=',wave,' silprof=',silprof
!***
      RETURN
      END


