      SUBROUTINE FE_ENTHALPY_STATE(A,UMIN,UMAX,NSTATE,U,UA,UB,METHOD,
     &                             T,DELTA_T,NSET,NMODES,EMODES,DE)

!********           fe_enthalpy_state v2       ********

      IMPLICIT NONE

! arguments

      CHARACTER METHOD*5
      INTEGER NSET,NMODES,NSTATE
      DOUBLE PRECISION A,DE,UMAX,UMIN
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
     &   HC,PI,TMINUS,TPLUS

      DOUBLE PRECISION EMODES_NEW(1000)

! --------------------------------------------------------------------
! subroutine FE_ENTHALPY_STATE
! given:
!     A          = radius (cm)
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
!     EMODES[]   = energies of vibrational modes (cm-1)

! returns:

!     U[]        = representative energy (erg) of bins 1-NSTATE
!     UA[],UB[]  = lower,upper limits (erg) to bins 1-NSTATE
!     T[]        = temperature (K) associated with each bin
!     DELTA_T[]  = temperature range (K) for each bin [*** NEVER USED ***]

! History:
! 11.02.19 (BTD) * adapted from sil_enthalpy_state
! 16.05.05 (BTD) v2
!                * modified calls to fe_temp [eliminate argument IMODE]
! end  history
!-----------------------------------------------------------------------
!*** diagnostic
!      write(0,*)'fe_enthalpy_state_v2 ckpt 0, method=',method
!      write(0,fmt='(a,1pe10.3)')'a=',a
!      write(0,fmt='(a,1pe10.3)')'umin=',umin
!      write(0,fmt='(a,1pe10.3)')'umax=',umax
!***
      PI=4.D0*ATAN(1.D0)
      HC=6.62607D-27*2.99792D10

      DENSITY=7.87
      NATOM=NINT((4.*PI/3.)*A**3*DENSITY/(55.845*1.6605D-24))

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
            IF((U(I)*(10.**DELTA_LGU)-U(I)).GT.(U(I)-U(I-1)))THEN
               JSET=I
               GOTO 555
            ENDIF
            IF(I.EQ.NSTATE)THEN
               WRITE(0,*)'fe_enthalpy_state: No Way!!!!'
               GOTO 555
            ENDIF
         ENDDO
 555     CONTINUE

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

!*** diagnostic
!      write(0,*)'fe_enthalpy_state_v2 ckpt 2'
!      write(0,*)'   u(1)=',u(1)
!      write(0,*)'   u(2)=',u(2)
!      write(0,*)'   u(3)=',u(3)
!      write(0,*)'   u(4)=',u(4)
!      write(0,*)'   u(5)=',u(5)
!***
! ----------------------------------------------------
! from u(i) --> T(i): for a given internal energy [cm-1], derive the 
! temperature of Fe cluster: EKT=E/(KT), 
! u(i): erg ---> energy (cm^{-1}) so divided by h*c;
! ck=C/K=dE/dT; T[i] == grain temperature; 
! modified on March 30, 2000: adopt BTD's subroutines.
! silicates: 3.5A-->natom=13.
! mode=-20 for silicates <= 25A; mode=1: Debye model;
! for simplicity mode can be set 1 in case of debye/dbcon.
! 08.03.22 (BTD) change to use MODE=0 for all sizes
!                when METHOD.NE.'stati'
      IF(METHOD.EQ.'stati')MODE=0
      IF(METHOD.NE.'stati')THEN
!         IF(A.LE.25.D-8)MODE=-20
!         IF(A.GT.25.D-8)MODE=1
         MODE=1
      ENDIF

!*** diagnostic
!      write(0,*)'fe_enthalpy_state_v2 ckpt 10, nstate=',nstate
!***
      DO I=1,NSTATE
         ENERGY=U(I)/HC
         IF(ENERGY.EQ.0.)THEN
            T(I)=0.
         ELSEIF(ENERGY.GT.0.)THEN

! 160505 (BTD) removed IMODE from FE_TEMP arg list

!*** diagnostic
!            write(0,*)'fe_enthalpy_state_v2 ckpt 11, energy=',energy
!***
            CALL FE_TEMP(ENERGY,NATOM,T(I),EKT,CK)
!*** diagnostic
!            write(0,fmt='(a,i3,1p2e10.3)')
!     &          'fe_enthalpy_state_v2 ckpt 12: i,energy,T(i)=',
!     &          I,ENERGY,T(I)
!***
         ENDIF

         ENERGYA=UA(I)/HC
         IF(ENERGYA.EQ.0.)THEN
            TMINUS=0.D0
         ELSEIF(ENERGYA.GT.0.)THEN

! 160505 (BTD) removed IMODE from FE_TEMP arg list

            CALL FE_TEMP(ENERGYA,NATOM,TMINUS,EKT,CK)
         ENDIF

         ENERGYB=UB(I)/HC
         IF(ENERGYB.EQ.0.)THEN
            TPLUS=0.D0
         ELSEIF(ENERGYB.GT.0.)THEN

! 160505 (BTD) removed IMODE from FE_TEMP arg list

            CALL FE_TEMP(ENERGYB,NATOM,TPLUS,EKT,CK)
         ENDIF

         DELTA_T(I)=TPLUS-TMINUS
	
      ENDDO

!*** diagnostic
!      write(0,*)'fe_enthalpy_state_v2 ckpt 9'
!      do i=1,nstate
!         write(0,fmt='(a,i3,1p2e10.3)')'i u(i)/hc, T(i)=',i,u(i)/hc,t(i)
!      enddo
!***
      RETURN
      END
