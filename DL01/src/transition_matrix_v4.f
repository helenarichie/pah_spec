      SUBROUTINE TRANSITION_MATRIX(METHOD,ICASE,BOVERA,AGR,
     &                             DTOH,NISRF,ISRF_WL,ISRF,DLGLAMBDA,
     &                             CABS,CSCA,NSTATE,U,UA,UB,LNGU,T,
     &                             AMATRIX,COOLING_TIME,RAD_TIME,
     &                             PSTATE,BGD)

      IMPLICIT NONE

!------------------------- transition_matrix_v4 -----------------------
! arguments

      CHARACTER METHOD*5
      INTEGER ICASE,NISRF,NSTATE
      DOUBLE PRECISION AGR,BOVERA,DLGLAMBDA,DTOH
      DOUBLE PRECISION
     &   AMATRIX(NSTATE,NSTATE),
     &   BGD(NSTATE,NSTATE),
     &   CABS(NISRF),
     &   COOLING_TIME(NSTATE),
     &   CSCA(NISRF),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF),
     &   LNGU(NSTATE),
     &   PSTATE(NSTATE),
     &   RAD_TIME(NSTATE),
     &   T(NSTATE),
     &   U(NSTATE),
     &   UA(NSTATE),
     &   UB(NSTATE)

! arguments of READQLIB

      REAL AEFF,WAVE
      REAL
     &   QABS(1:3),
     &   QEXT(1:3),
     &   QSCA(1:3)

! local variables

      INTEGER I,F,J,K
      REAL BOVERA_SP
      DOUBLE PRECISION
     &   ADIAG,AFI_HEATING,AFI_COOLING,AMATRIX_DBCON,
     &   COOLING_RATE,PI,RAD_RATE,SUM_PSTATE

!-----------------------------------------------------------------------
! subroutine TRANSITION_MATRIX v4
! given
!     METHOD = 'stati' for exact-statistical treatment
!              'dbdis' for Debye-discrete-thermal treatment
!              'dbcon' for Debye-continuous-thermal treatment
!              'dbdgl' for ???
!     AGR     = grain radius (cm) [to calculate C_abs]
!     DTOH    = D/H ratio (used only by PAH option)
!     NISRF   = number of wavelengths in tabulated ISRF
!     ISRF_WL[]= wavelengths (cm) for tabulated ISRF
!                (assumed to be uniform in log(lambda))
!     ISRF[]  = c*u_lambda (erg cm-3 s-1) for tabulated ISRF
!     DLGLAMBDA = delta(log10(lambda)) for ISRF_WL
!     CABS[]  = C_abs (cm2) at wavelengths ISRF_WL
!     CSCA[] = C_sca (cm2) at wavelengths ISRF_WL
!     NSTATE  = number of energy bins to use for transition matrix
!     U[J]    = mid-point of bin J (erg)
!     UA[J]   = lower limit of bin J (erg)
!     UB[J]   = upper limit of bin J (erg)
!     LNGU[J] = ln(g) for bin J, g=degeneracy
!     T[J]    = "temperature" associated with bin J

! returns

!     AMATRIX[J,I] = transition rate (s-1) from bin I to bin J
!     COOLING_TIME[J]=U/(dU/dt)_cool (sec) for bin J
!     RAD_TIME[J]  = 1/(photon absorption rate) for bin J
!     PSTATE[J] = preliminary estimate for steady state probability 
!                 of being in bin J .
!                 This is NOT the final result, but is produced only for
!                 preconditioning for BiCGstab solution method.

!-----------------------------------------------------------------------
! Subroutine TRANSITION_MATRIX calculates the transition matrix
! A_fi, including contributions from absorption of photons, and
! emission of photons.  It also produces an approximate estimate
! (using method of Guhathakurta & Draine) for the probability
! distribution PSTATE.

! Routine originally written by Aigen Li, Princeton University
! History:
! 99.09.28 (AL)  original version
! 99.10.18 (AL)  modifications
! 99.10.28 (AL)  include silicates
! 00.01.27 (AL)  modified
! 00.02.18 (AL)  include statistical and thermal treatments
! 00.03.08 (AL)  add DTOH = D/H for PAHs (DTOH=0 for sil/gra)
! 00.11.12 (BTD) cosmetic changes, added comments
! 07.12.10 (BTD) changed type of ZG to INTEGER
! 13.07.12 (BTD) v2: replace calls to SIL_GRA_PAH_CRS_SCT
!                    with call to new routine READQLIB
! 14.02.03 (BTD) typed and defined PIA2
! 14.02.17 (BTD) v3
!                * add ICASE,ISHAP to argument list
!                * add ISHAP to argument list of READQLIB
!                * remove local calculation of ICASE since it
!                  is now being input as argument
!                * eliminate variable COMPOSITION
!                * eliminate variable ZG
! 17.10.10 (BTD) * eliminate ISHAP
!                * introduce BOVERA to argument list
! 18.08.13 (BTD) v4
!                * remove radius from many arg lists where not needed
! end history
!-----------------------------------------------------------------------
! calculate the absorption cross section;
! BE CAREFUL: here isrf_wl(i) -- cm;
! AGR: cm; Cabs: cm^2;

!*** diagnostic
      write(0,fmt='(a,i3)')'transition_matrix_v4 ckpt 1: icase=',icase
!      write(0,*)'  i    ua       u       ub'
!      do i=1,nstate
!         write(0,fmt='(i5,1p3e10.3)')i,ua(i),u(i),ub(i)
!      enddo
!***

!----------------------------------------------------------------------
! Prepare for calls to READQLIB:
! convert from double precision and cm to single precision and um
      PI=4.D0*ATAN(1.D0)

      AEFF=1.E4*REAL(AGR)   ! aeff (um)
      BOVERA_SP=REAL(BOVERA)
!-----------------------------------------------------------------------
      DO I=1,NISRF
         WAVE=1.E4*REAL(ISRF_WL(I))
         CALL READQLIB(ICASE,BOVERA_SP,AEFF,WAVE,QABS,QEXT,QSCA)
         CABS(I)=PI*AGR**2*(QABS(1)+QABS(2)+QABS(3))/3.
      ENDDO

!*** diagnostic
!      write(0,*)'transition_matrix_v4 ckpt 2'
!      write(0,*)'  returned from sil_gra_pah_crs_sct'
!      write(0,*)' begin looping over I and F ...'
!***
! calculate Afi

      DO I=1,NSTATE
         DO F=1,NSTATE
!*** diagnostic
!      write(0,fmt='(a,i4,i4)')'transition_matrix_v4 ckpt 3: I,F=',I,F
!***
            AMATRIX(F,I)=0.

            IF(I.LT.F)THEN
!*** diagnostic
!      write(0,*)'transition_matrix_v4 ckpt 4'
!      write(0,*)' call heating_afi ...'
!***
               CALL HEATING_AFI(I,U(I),UA(I),UB(I),F,U(F),UA(F),UB(F),
     &                          NSTATE,NISRF,ISRF_WL,ISRF,CABS,
     &                          AFI_HEATING)
!*** diagnostic
!      write(0,*)'transition_matrix_v4 ckpt 5'
!      write(0,*)' returned from heating_afi'
!***
               AMATRIX(F,I)=AFI_HEATING

            ELSEIF(I.GT.F)THEN
               IF(METHOD.EQ.'stati')THEN
!*** diagnostic
!      write(0,*)'transition_matrix_v4 ckpt 6'
!      write(0,*)' call stat_cooling_afi ...'
!***
                  CALL STAT_COOLING_AFI(I,U(I),UA(I),UB(I),LNGU(I),
     &                                  F,U(F),UA(F),UB(F),LNGU(F),
     &                                  NSTATE,NISRF,ISRF_WL,ISRF,
     &                                  CABS,AFI_COOLING)
!*** diagnostic
!      write(0,*)'transition_matrix_v4 ckpt 7'
!      write(0,*)' returned from stat_cooling_afi'
!***
               ELSEIF(METHOD.EQ.'dbdis'.OR.METHOD.EQ.'dbcon')THEN
!*** diagnostic
!                 write(0,*)'transition_matrix_v4 ckpt 8'
!***
                 CALL THERM_COOLING_AFI(I,T(I),U(I),UA(I),UB(I),
     &                                   F,U(F),UA(F),UB(F),
     &                                   (UB(2)-UA(2)),NSTATE,
     &                                   NISRF,ISRF_WL,ISRF,CABS,
     &                                   AFI_COOLING)
!*** diagnostic
!                 write(0,*)'transition_matrix_v4 ckpt 9'
!***

               ELSEIF(METHOD.EQ.'dbgdl')THEN
                  IF(I.EQ.2)THEN
                     CALL THERM_COOLING_AFI(I,T(I),U(I),UA(I),UB(I),
     &                                      F,U(F),UA(F),UB(F),
     &                                      (UB(2)-UA(2)),NSTATE,
     &                                      NISRF,ISRF_WL,ISRF,CABS,
     &                                      AFI_COOLING)
                  ELSEIF(I.NE.2)THEN
                     CALL DBGDL_COOLING_AFI(I,T(I),U(I),UA(I),UB(I),
     &                                      F,U(F),UA(F),UB(F),
     &                                      NSTATE,
     &                                      NISRF,ISRF_WL,ISRF,CABS,
     &                                      AFI_COOLING)
                  ENDIF  
               ENDIF
               AMATRIX(F,I)=AFI_COOLING

            ELSEIF(I.EQ.F)THEN
               AMATRIX(F,I)=0.D0
            ENDIF

         ENDDO
      ENDDO

!*** diagnostic
!      write(0,*)'transition_matrix_v4 ckpt 10'
!      write(0,*)' completed loops over I and F'
!***

! the following lines are only for the continuous cooling method!
! when call dbgdl_cooling_Afi for the continuous cooling method,
! should comment out the following lines. if call therm_cooling_Afi
! for the continuous cooling method, the following lines are 
! required to construct A[i-1,i]. March 03, 2000.

      IF(METHOD.EQ.'dbcon')THEN
         DO I=2,NSTATE
            AMATRIX_DBCON=0.D0
            DO F=I-1,1,-1
               AMATRIX_DBCON=AMATRIX_DBCON+AMATRIX(F,I)*(U(I)-U(F))
               IF(I.NE.(F+1))AMATRIX(F,I)=0.
            ENDDO
            AMATRIX(I-1,I)=AMATRIX_DBCON/(U(I)-U(I-1))
         ENDDO
      ENDIF

! the above lines are only for the continuous cooling method!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! calculate Aii:
! ==============
      DO I=1,NSTATE 
         ADIAG=0.D0
         DO F=1,NSTATE
            ADIAG=ADIAG+AMATRIX(F,I)
         ENDDO
         AMATRIX(I,I)=-ADIAG
      ENDDO

!        WRITE(0,*)'The transition matrix has been constructed!'
!       =======================================================

! -------------------------------------------------------------
! Estimate the cooling time:

      DO I=1,NSTATE		
         IF(I.EQ.1)THEN
            COOLING_TIME(I)=0.
            RAD_TIME(I)=0.
         ELSEIF(I.NE.1)THEN

	    COOLING_RATE=0.
            RAD_RATE=0.
            DO K=1,I-1
               COOLING_RATE=COOLING_RATE+AMATRIX(K,I)*(U(I)-U(K))
               RAD_RATE=RAD_RATE+AMATRIX(K,I)
            ENDDO

            COOLING_TIME(I)=U(I)/COOLING_RATE
            RAD_TIME(I)=1.D0/RAD_RATE
         ENDIF	
      ENDDO	

! ------------------------------------------------------------

      IF(METHOD.NE.'dbcon' .AND. METHOD.NE.'dbgdl')THEN
         DO F=1,NSTATE
            PSTATE(F)=0.D0
            DO J=1,NSTATE
               BGD(F,J)=0.D0
            ENDDO
         ENDDO
         GOTO 100
      ELSEIF(METHOD.EQ.'DBCON' .OR. METHOD.EQ.'DBGDL')THEN
         GOTO 200
      ENDIF

 200  CONTINUE

! Adopt the method of Guhathakurta & Draine (1989)
! to estimate the temperature probability distribution;
! the results will be used as pre-conditioning input in
! calling BiCGstab subroutine.
 
      PSTATE(1)=1.
      DO F=2,NSTATE
         PSTATE(F)=0.
         DO J=1,F-1
            BGD(F,J)=0.
            DO K=F,NSTATE
               BGD(F,J)=BGD(F,J)+AMATRIX(K,J)
            ENDDO
            PSTATE(F)=PSTATE(F)+BGD(F,J)/AMATRIX(F-1,F)*PSTATE(J)
         ENDDO
      ENDDO

! normalize to 1:

      SUM_PSTATE=0.
      DO I=1,NSTATE
         SUM_PSTATE=SUM_PSTATE+PSTATE(I)
      ENDDO

      DO I=1,NSTATE
         PSTATE(I)=PSTATE(I)/SUM_PSTATE
      ENDDO

 100  CONTINUE
!*** diagnostic
!      write(0,*)'transition_matrix_v4 ckpt 99'
!***
      RETURN
      END
