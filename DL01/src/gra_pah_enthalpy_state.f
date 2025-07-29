      SUBROUTINE GRA_PAH_ENTHALPY_STATE(NC,NH,ND,UMIN,UMAX,NSTATE,
     &                                  U,UA,UB,METHOD,T,DELTA_T,
     &                                  NSET,NMODES,EMODES,DE)

      IMPLICIT NONE

c arguments:

      CHARACTER METHOD*5	
      INTEGER NC,ND,NH,NSTATE,NSET,NMODES
      DOUBLE PRECISION DE,UMIN,UMAX
      DOUBLE PRECISION
     &   DELTA_T(NSTATE),
     &   EMODES(NMODES),
     &   T(NSTATE),
     &   U(NSTATE),
     &   UA(NSTATE),
     &   UB(NSTATE)

c local variables:

      INTEGER I,JSET,MODE,NSET_MAX
      DOUBLE PRECISION 
     &   CK,DELTA_LGU,EKT,ENERGY,ENERGYA,ENERGYB,HC,TMINUS,TPLUS

      DOUBLE PRECISION EMODES_NEW(1000)
c-----------------------------------------------------------------------
c subroutine GRA_PAH_ENTHALPY_STATE
c given:
c     METHOD   = 'stati' to ???
c     NC       = number of C atoms
c     ND       = number of D atoms
c     NH       = number of H atoms (not including D)
c     NSTATE   = number of energy states desired
c     NMODES   = number of vibrational modes
c     EMODES[] = energies of vibrational modes (cm-1)
c     UMIN     = lowest energy of interest (erg)
c     UMAX     = highest energy of interest (erg)
c     DE       = energy increment (cm-1) used in Beyer-Swinehart
c                algorithm calculation of number of vibrational
c                states.  We require that enthalpy bin boundaries
c                be integral multiples of DE to avoid interpolation
c                problems when calculating degeneracy of each bin.
c   
c returns
c     NSET     =
c     U[J]     = mid-bin energy for bin J (erg)
c     UA[J]    = lower bound energy for bin J (erg)
c     UB[J]    = upper bound energy for bin J (erg)
c
c-----------------------------------------------------------------------
c Basic purpose is to create a range of energy bins for PAHs
c Original version written by Aigen Li, Princeton University
c Modified by B.T. Draine, Princeton University
c History
c 00.01.27 (AL)
c 00.02.21 (AL)  modified
c 00.02.28 (AL)  modified
c 00.03.08 (AL)  modified to allow deuterium
c 00.11.28 (BTD) reduced EMODES_NEW to array of length 1000
c 06.05.23 (BTD) changed to use MODE=0 for method.ne.stati and
c                NC.le.7360 
c                (previously was using mode=-20)
c 10.03.06 (BTD) modified setting of bins to ensure that
c                bin boundaries are monotonic,
c                nominal bin energy is now set to average of
c                bin boundaries
c end history
c-----------------------------------------------------------------------

      HC=6.62607D-27*2.99792D10

c*** diagnostic
c      write(0,*)'gra_pah_enthalpy_state ckpt 0'
c      write(0,*)'    METHOD=',method
c      write(0,*)'     UMIN =',umin
c      write(0,*)'       DE =',DE
c      write(0,*)' EMODES(1)=',EMODES(1)
c      write(0,*)' EMODES(2)=',EMODES(2)
c      write(0,*)' EMODES(3)=',EMODES(3)
c      write(0,*)' EMODES(4)=',EMODES(4)
c      write(0,*)' EMODES(5)=',EMODES(5)
c      write(0,*)' EMODES(6)=',EMODES(6)
c      write(0,*)' EMODES(7)=',EMODES(7)
c***

      IF(UMIN.EQ.0.)THEN 

c bin 1 = ground state (zero width)

         U(1)=0.D0
         UA(1)=0.D0
         UB(1)=0.D0

c in case of DBGDL (continuous cooling Debye-thermal model),
c we do not need to construct energy levels specifically, but
c just logarithmly equal-space them. 00.03.07.
c simplified on April 04, 2000.
c modified on July 11, 2000: take u[1]=0, u[2]=E1, u[3]=E2, 
c u[4]=E4, u[5]=E6, u[i]=E(i*2-4).... u[nset];

c construct EMODES_NEW(1,2,...,NSET_MAX) = distinct mode energies
c in order of increasing frequency.
c EMODES_NEW(1) = EMODES(1) = first vibrational state

         NSET_MAX=0
         DO I=1,NMODES-1
            IF(EMODES(I+1).NE.EMODES(I))THEN
               NSET_MAX=NSET_MAX+1
               EMODES_NEW(NSET_MAX)=EMODES(I)
               IF(NSET_MAX.EQ.1000)GOTO 1000
            ENDIF
         ENDDO

c*** diagnostic
c         write(0,*)'gra_pah_enthalpy_state ckpt 10'
c         write(0,*)'  emodes_new(1)=',emodes_new(1)
c         write(0,*)'  emodes_new(2)=',emodes_new(2)
c         write(0,*)'  emodes_new(3)=',emodes_new(3)
c         write(0,*)'  emodes_new(4)=',emodes_new(4)
c         write(0,*)'  emodes_new(5)=',emodes_new(5)
c         write(0,*)'  emodes_new(6)=',emodes_new(6)
c         write(0,*)'  emodes_new(7)=',emodes_new(7)
c***
 1000    NSET=(NSET_MAX+4)/2
         IF(NSET.GE.11)NSET=11

c note nset*2-4 <= nset_max! (for C24H12 nset<=36!; 
c for a=3.5A, nc=20, nh=10, nset<=30);
c on July 19, 2000: set nset=11 (=18 modes).
c note Bruce's 10 implies 11 here!

c bin 2:
         UA(2)=(3.*EMODES_NEW(1)-EMODES_NEW(2))/2.
         UB(2)=(EMODES_NEW(1)+EMODES_NEW(2))/2.
         U(2)=(UA(2)+UB(2))/2.

c bin 3:
         UA(3)=UB(2)
         UB(3)=(EMODES_NEW(2)+EMODES_NEW(3))/2.	
         U(3)=(UA(3)+UB(3))/2.

c bins 4 - NSET each contain two normal modes (each of which
c               may possibly be degenerate):

         DO I=4,NSET
            UB(I)=(EMODES_NEW(I*2-4)+EMODES_NEW(I*2-3))/2.	
            UA(I)=UB(I-1)
            U(I)=(UA(I)+UB(I))/2.
         ENDDO

c*** diagnostic
c         write(0,*)'gra_pah_enthalpy_state ckpt 15'
c         write(0,*)'  nset=',nset
c         write(0,*)' ua(1)=',ua(1)
c         write(0,*)'  u(1)=',u(1)
c         write(0,*)' ub(1)=',ub(1)
c         write(0,*)' ua(2)=',ua(2)
c         write(0,*)'  u(2)=',u(2)
c         write(0,*)' ub(2)=',ub(2)
c         write(0,*)' ua(3)=',ua(3)
c         write(0,*)'  u(3)=',u(3)
c         write(0,*)' ub(3)=',ub(3)
c         write(0,*)' ua(4)=',ua(4)
c         write(0,*)'  u(4)=',u(4)
c         write(0,*)' ub(4)=',ub(4)
c         write(0,*)' ua(5)=',ua(5)
c         write(0,*)'  u(5)=',u(5)
c         write(0,*)' ub(5)=',ub(5)
c***

c beyond NSET, bins are first initially uniformly space in E,
c              and then switch over to uniform spacing in log(E)

         DO I=NSET+1,NSTATE
            U(I)=U(NSET)+(UB(NSET)-UA(NSET))*(I-NSET)
            DELTA_LGU=DLOG10((UMAX/HC)/U(I))/(NSTATE-I)
            IF((U(I)*(10.**DELTA_LGU)-U(I)).GT.(U(I)-U(I-1)))THEN
               JSET=I
               GOTO 555
            ENDIF
            IF(I.EQ.NSTATE)THEN
               WRITE(0,*)'No Way -- this cannot be right!'
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

c*** diagnostic
c         write(0,*)'gra_pah_enthalpy_state ckpt 16'
c         write(0,*)'  nset=',nset
c         write(0,*)' ua(1)=',ua(1)
c         write(0,*)'  u(1)=',u(1)
c         write(0,*)' ub(1)=',ub(1)
c         write(0,*)' ua(2)=',ua(2)
c         write(0,*)'  u(2)=',u(2)
c         write(0,*)' ub(2)=',ub(2)
c         write(0,*)' ua(3)=',ua(3)
c         write(0,*)'  u(3)=',u(3)
c         write(0,*)' ub(3)=',ub(3)
c         write(0,*)' ua(4)=',ua(4)
c         write(0,*)'  u(4)=',u(4)
c         write(0,*)' ub(4)=',ub(4)
c         write(0,*)' ua(5)=',ua(5)
c         write(0,*)'  u(5)=',u(5)
c         write(0,*)' ub(5)=',ub(5)
c         write(0,*)' ua(6)=',ua(6)
c         write(0,*)'  u(6)=',u(6)
c         write(0,*)' ub(6)=',ub(6)
c         write(0,*)' ua(7)=',ua(7)
c         write(0,*)'  u(7)=',u(7)
c         write(0,*)' ub(7)=',ub(7)
c***

c now set bin boundaries at multiple of DE
c but make sure that bin boundaries are monotonic

         DO I=1,NSTATE
            UA(I)=DE*NINT(UA(I)/DE)
            UB(I)=DE*NINT(UB(I)/DE)
         ENDDO
         DO I=2,NSTATE
            UA(I)=UB(I-1)
            IF(UB(I).LT.UA(I)+DE)UB(I)=UA(I)+DE
         ENDDO
         DO I=2,NSTATE
            U(I)=0.5*(UA(I)+UB(I))
         ENDDO

c*** diagnostic
c         write(0,*)'gra_pah_enthalpy_state ckpt 17'
c         write(0,*)'  nset=',nset
c         write(0,*)' ua(1)=',ua(1)
c         write(0,*)'  u(1)=',u(1)
c         write(0,*)' ub(1)=',ub(1)
c         write(0,*)' ua(2)=',ua(2)
c         write(0,*)'  u(2)=',u(2)
c         write(0,*)' ub(2)=',ub(2)
c         write(0,*)' ua(3)=',ua(3)
c         write(0,*)'  u(3)=',u(3)
c         write(0,*)' ub(3)=',ub(3)
c         write(0,*)' ua(4)=',ua(4)
c         write(0,*)'  u(4)=',u(4)
c         write(0,*)' ub(4)=',ub(4)
c         write(0,*)' ua(5)=',ua(5)
c         write(0,*)'  u(5)=',u(5)
c         write(0,*)' ub(5)=',ub(5)
c         write(0,*)' ua(6)=',ua(6)
c         write(0,*)'  u(6)=',u(6)
c         write(0,*)' ub(6)=',ub(6)
c         write(0,*)' ua(7)=',ua(7)
c         write(0,*)'  u(7)=',u(7)
c         write(0,*)' ub(7)=',ub(7)
c***

c convert U,UA,UB from cm-1 to erg

         DO I=1,NSTATE
            U(I)=U(I)*HC
            UA(I)=UA(I)*HC
            UB(I)=UB(I)*HC
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

c note umin,umax (thus u,ua,ub) are already in erg!

      ENDIF

c ----------------------------------------------------		
c from u(i) --> T(i): for a given enthalpy, derive the 
c temperature of PAHs; energy: cm^{-1}; ekt=E/(KT), 
c u(i): erg ---> energy (cm^{-1}) so devided by h*c;
c ck=C/K=dE/dT; T[i] == grain temperature; 
c modified on Feb.18, 2000: add ``mode'' as argument.
c mode=-20 for a<=30A; mode=1 for C-C Debye model.
c for simplicity can set mode=1 in case of debye/dbcon.
c nc=7361 ->a=20A for PAH (2.24g/cm^3).
 
      IF(METHOD.EQ.'STATI')MODE=0
      IF(METHOD.NE.'STATI')THEN

c-----------------------------------------------------------------------
c 06.05.23 (BTD) change from using MODE=-20 to using MODE=0
c	 IF(NC.LE.7360)MODE=-20
	 IF(NC.LE.7360)MODE=0
c-----------------------------------------------------------------------
	 IF(NC.GT.7360)MODE=1
      ENDIF
c*** diagnostic
c      write(0,*)'gra_pah_enthalpy_state ckpt 20'
c      write(0,*)' ua(1)=',ua(1)
c      write(0,*)'  u(1)=',u(1)
c      write(0,*)' ub(1)=',ub(1)
c      write(0,*)' ua(2)=',ua(2)
c      write(0,*)'  u(2)=',u(2)
c      write(0,*)' ub(2)=',ub(2)
c      write(0,*)' ua(3)=',ua(3)
c      write(0,*)'  u(3)=',u(3)
c      write(0,*)' ub(3)=',ub(3)
c      write(0,*)' ua(4)=',ua(4)
c      write(0,*)'  u(4)=',u(4)
c      write(0,*)' ub(4)=',ub(4)
c      write(0,*)' ua(5)=',ua(5)
c      write(0,*)'  u(5)=',u(5)
c      write(0,*)' ub(5)=',ub(5)
c***

      DO I=1,NSTATE
         ENERGY=U(I)/HC
         IF(ENERGY.EQ.0.)THEN
            T(I)=0.
         ELSEIF(ENERGY.GT.0.)THEN
            CALL PAH_TEMP(MODE,ENERGY,NC,NH,ND,T(I),EKT,CK)
         ENDIF

         ENERGYA=UA(I)/HC
         IF(ENERGYA.EQ.0.)THEN
            TMINUS=0.
         ELSEIF(ENERGYA.GT.0.)THEN
            CALL PAH_TEMP(MODE,ENERGYA,NC,NH,ND,TMINUS,EKT,CK)	
         ENDIF

         ENERGYB=UB(I)/HC
         IF(ENERGYB.EQ.0.)THEN
            TPLUS=0.
         ELSEIF(ENERGYB.GT.0.)THEN
            CALL PAH_TEMP(MODE,ENERGYB,NC,NH,ND,TPLUS,EKT,CK)
         ENDIF

         DELTA_T(I)=TPLUS-TMINUS

c*** diagnostic
c    sanity check
         if(i.gt.1)then
            if(t(i).le.t(i-1))then
               write(0,*)'gra_pah_enthalpy_state ckpt 98'
               write(0,*)'   pah_temp called with mode=',mode
               write(0,*)'   energya=',energya
               write(0,*)'   energyb=',energyb
               write(0,*)'   fatal error in gra_pah_enthalpy_state:'
               write(0,*)'   i=',i,' t(i-1)=',t(i-1),' > t(i)=',t(i)
               stop
            endif
            if(delta_t(i).le.0.)then
               write(0,*)'gra_pah_enthalpy_state ckpt 99'
               write(0,*)'   pah_temp called with mode=',mode
               write(0,*)'   energya=',energya,' tminus=',tminus
               write(0,*)'   energy =',energy,'   t(i)=',t(i)
               write(0,*)'   energyb=',energyb,'  tplus=',tplus
               write(0,*)'   fatal error in gra_pah_enthalpy_state:'
               write(0,*)'   i=',i,' delta_t(i)=',delta_t(i),' < 0'
               stop
            endif
         endif
c***
      ENDDO

      RETURN
      END
