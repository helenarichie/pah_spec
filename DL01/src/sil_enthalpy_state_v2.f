      SUBROUTINE SIL_ENTHALPY_STATE(AESVS,UMIN,UMAX,NSTATE,U,UA,UB,
     &                              METHOD,T,DELTA_T,NSET,NMODES,
     &                              EMODES,DE)

      IMPLICIT NONE
!                       sil_enthalpy_state_v2
! arguments

      CHARACTER METHOD*5
      INTEGER NSET,NMODES,NSTATE
      DOUBLE PRECISION AESVS,DE,UMAX,UMIN
      DOUBLE PRECISION
     &   DELTA_T(NSTATE),
     &   EMODES(NMODES),
     &   T(NSTATE),
     &   U(NSTATE),
     &   UA(NSTATE),
     &   UB(NSTATE)

! local variables

      INTEGER I,JSET,MODE,NATOM,NSET_MAX
      DOUBLE PRECISION 
     &   CK,DELTA_LGU,DENSITY,EKT,ENERGY,ENERGYA,ENERGYB,
     &   HC,MAW,TMINUS,TPLUS

      DOUBLE PRECISION EMODES_NEW(1000)

!-----------------------------------------------------------------------
! subroutine SIL_ENTHALPY_STATE
! given:
!     AESVS      = equal-solid-volume-sphere radius (cm)
!                = a_eff*(1-P)^{1/3} where P = porosity
!     UMIN       = lowest energy of interest (erg)
!     UMAX       = highest energy of interest (erg)
!     NSTATE     = number of energy bins desired
!     METHOD     = 'stati' for exact-statistical method
!                  'dbdis' or 'dbcon' for thermal methods
!     NMODES     = number of vibrational modes
!     DE         = energy increment (cm-1) used in Beyer-Swinehart method
!                  We require that enthalpy bin boundaries be integral
!                  multiples of DE to avoid interpolation problems when
!                  calculating degeneracy of each bin.
!     EMODES[J]  = energies of vibrational modes (cm-1)
!                  J=1-NMODES

! returns:

!     U[]        = representative energy (erg) of bins 1-NSTATE
!     UA[],UB[]  = lower,upper limits (erg) to bins 1-NSTATE
!     T[]        = temperature (K) associated with each bin
!     DELTA_T[]  = temperature range (K) for each bin [*** NEVER USED ***]
!     NSET       =

! Originally written by Aigen Li, Princeton University
! History:
! 99.10.21 (AL)  v1 First written
! 00.03.20 (AL)  Modified after Taiwan's election.
! 00.11.26 (BTD) Cosmetic changes; added comments
! 00.11.28 (BTD) reduced length of EMODES_NEW -- 1000 is more than
!                sufficient
! 08.03.22 (BTD) Changed to use MODE=0 for all sizes when determining
!                "Temperature" of a level.
! 18.08.11 (BTD) v2 
!                * modified to clarify treatment of porous grains
!                * now input AESVS = radius of equal-mass sphere
! 22.12.07 (BTD) * Cosmetics
!                * change silicate mass density to 3.41 g cm-3
!                  as estimated in Draine & Hensley 2021a
!                  "The Dielectric Function of Astrodust..." ApJ 909:94
! end  history
!-----------------------------------------------------------------------

      HC=6.62607D-27*2.99792D10

!      DENSITY=3.8    ! previous estimate for silicate mass density

      DENSITY=3.41

! MAW = mean atomic weight
!     = 172./7.   = 24.6 amu for MgFeSiO4
!     = 133.6/6.2 = 21.5 amu for Mg_1.3 Fe_0.3 Si O_3.6

      MAW=21.5
      NATOM=NINT((4.*3.1416/3.)*AESVS**3*DENSITY/(MAW*1.66D-24))

!*** diagnostic
      write(0,fmt='(a,i10)')'sil_enthalpy_state ckpt 1: natom=',natom
      write(0,fmt='(a,i10)')'     nset=',nset
      write(0,fmt='(a,i10)')'   nmodes=',nmodes
      write(0,fmt='(a,1pe10.3)')'     emodes(1)=',emodes(1)
      write(0,fmt='(a,1pe10.3)')'emodes(nmodes)=',emodes(nmodes)
!***
      IF(UMIN.EQ.0.)THEN

! bin 1 = ground state (zero width)

         U(1)=0.D0
         UA(1)=0.D0
         UB(1)=0.D0

         NSET_MAX=0
         DO I=1,NMODES-1
            IF(EMODES(I+1).NE.EMODES(I))THEN
               NSET_MAX=NSET_MAX+1
               EMODES_NEW(NSET_MAX)=EMODES(I)
               IF(NSET_MAX.EQ.1000)GOTO 1000
            ENDIF
         ENDDO

 1000    NSET=(NSET_MAX+4)/2
         IF(NSET.GE.11)NSET=11

!*** diagnostic
         write(0,fmt='(a,i8)')'sil_enthalpy_state ckpt 2, nset=',nset
!***
! bin 2:
         UA(2)=(3.*EMODES_NEW(1)-EMODES_NEW(2))/2.
         UB(2)=(EMODES_NEW(1)+EMODES_NEW(2))/2.
         U(2)=(UA(2)+UB(2))/2.

! bin 3:
         UA(3)=UB(2)
         UB(3)=(EMODES_NEW(2)+EMODES_NEW(3))/2.	
         U(3)=(UA(3)+UB(3))/2.

! bins 4 - NSET each contain two normal modes (each of which may
!               possibly be degenerate):

         DO I=4,NSET
            UB(I)=(EMODES_NEW(I*2-4)+EMODES_NEW(I*2-3))/2.	
            UA(I)=UB(I-1)
            U(I)=(UA(I)+UB(I))/2.
         ENDDO

! beyond NSET, bins are first uniformly spaced in E,
!              and then switch over to uniform spacing in log(E)

         DO I=NSET+1,NSTATE
            U(I)=U(NSET)+(UB(NSET)-UA(NSET))*(I-NSET)
            DELTA_LGU=DLOG10((UMAX/HC)/U(I))/(NSTATE-I)
!*** diagnostic
!            write(0,fmt='(a,i10)')
!     &         'sil_enthalpy_state ckpt 3: I=',I
!            write(0,fmt='(a,1pe10.3)')'   Delta_lgU=',DELTA_LGU
!***
            IF((U(I)*(10.**DELTA_LGU)-U(I)).GT.(U(I)-U(I-1)))THEN
               JSET=I
               GOTO 555
            ENDIF
            IF(I.EQ.NSTATE)THEN
               WRITE(0,*)'No Way!!!!'
               WRITE(0,*)'sil_enthalpy_state ckpt 4 fatal error'
               STOP
!               GOTO 555
            ENDIF
         ENDDO
 555     CONTINUE

!*** diagnostic
!         write(0,fmt='(a,i10)')'sil_enthalpy_state ckpt 5, NSET=',NSET
!         write(0,fmt='(a,1p2e10.3)')
!     &      '   ua(nset),ub(nset)=',ua(nset),ub(nset)
!***

         DO I=NSET+1,JSET
            U(I)=U(NSET)+(UB(NSET)-UA(NSET))*(I-NSET)
         ENDDO

         DO I=JSET+1,NSTATE
            U(I)=U(JSET)*10.**((I-JSET)*DELTA_LGU)
         ENDDO

         DO I=NSET+1,NSTATE-1
            UA(I)=(U(I)+U(I-1))/2.
            UB(I)=(U(I)+U(I+1))/2.
         ENDDO

         UA(NSTATE)=UB(NSTATE-1)
         UB(NSTATE)=U(NSTATE)

! convert U,UA,UB from cm-1 to erg

         DO I=1,NSTATE
            U(I)=U(I)*HC
            UA(I)=DE*NINT(UA(I)/DE)*HC
            UB(I)=DE*NINT(UB(I)/DE)*HC
         ENDDO

      ELSEIF(UMIN.NE.0.)THEN
         DELTA_LGU=DLOG10(UMAX/UMIN)/(NSTATE-1.)
!*** diagnostic
!         write(0,fmt='(a,1pe10.3)')
!     &      'sil_enthalpy_state ckpt 6, Delta_lgU=',DELTA_LGU
!***
         NSET=0
         DO I=1,NSTATE
            U(I)=UMIN*10.**((I-1.)*DELTA_LGU)
         ENDDO
         DO I=2,NSTATE-1	
            UA(I)=(U(I)+U(I-1))/2.
            UB(I)=(U(I)+U(I+1))/2.
         ENDDO
         UA(1)=U(1)
         UB(1)=UA(2)
         UA(NSTATE)=UB(NSTATE-1)
         UB(NSTATE)=U(NSTATE)

! note umin,umax (thus u,ua,ub) are already in erg!

      ENDIF    

! ----------------------------------------------------
! from u(i) --> T(i): for a given enthalpy, derive the 
! temperature of silicates; energy: cm^{-1}; ekt=E/(KT), 
! u(i): erg ---> energy (cm^{-1}) so devided by h*c;
! ck=C/K=dE/dT; T[i] == grain temperature; 
! modified on March 30, 2000: adopt BTD's subroutines.
! silicates: 3.5A-->natom=13.
! mode=-20 for silicates <= 25A; mode=1: Debye model;
! for simplicity mode can be set 1 in case of debye/dbcon.
! 08.03.22 (BTD) change to use MODE=0 for all sizes
!                when METHOD.NE.'stati'
      IF(METHOD.EQ.'stati')MODE=0
      IF(METHOD.NE.'stati')THEN
!         IF(AESVS.LE.25.D-8)MODE=-20
!         IF(AESVS.GT.25.D-8)MODE=1
         MODE=1
      ENDIF

      DO I=1,NSTATE
         ENERGY=U(I)/HC
         IF(ENERGY.EQ.0.)THEN
            T(I)=0.
         ELSEIF(ENERGY.GT.0.)THEN
            CALL SIL_TEMP(MODE,ENERGY,NATOM,T(I),EKT,CK)
         ENDIF

         ENERGYA=UA(I)/HC
         IF(ENERGYA.EQ.0.)THEN
            TMINUS=0.D0
         ELSEIF(ENERGYA.GT.0.)THEN
            CALL SIL_TEMP(MODE,ENERGYA,NATOM,TMINUS,EKT,CK)
         ENDIF

         ENERGYB=UB(I)/HC
         IF(ENERGYB.EQ.0.)THEN
            TPLUS=0.D0
         ELSEIF(ENERGYB.GT.0.)THEN
            CALL SIL_TEMP(MODE,ENERGYB,NATOM,TPLUS,EKT,CK)
         ENDIF

         DELTA_T(I)=TPLUS-TMINUS
	
      ENDDO

      RETURN
      END
