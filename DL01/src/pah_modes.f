      SUBROUTINE PAH_MODES(MXMODES,NC,NH,ND,NMODES,EMODES)
      IMPLICIT NONE
!
! Arguments:
!
      INTEGER MXMODES
      INTEGER NC,ND,NH,NMODES
      DOUBLE PRECISION EMODES(1:MXMODES)
!
! Common:
!
      DOUBLE PRECISION EV1_COM
      COMMON/EVCOM/EV1_COM
!
! Local variables:
!
      INTEGER J,N1,N_MD,N_S
      DOUBLE PRECISION BETA,DJOFF,EMCC_IP,EMCC_OP,EMCH_IP,EMCH_OP,
     &                 EMCH_ST,EMDH_IP,EMDH_OP,EMDH_ST
      SAVE EMCC_IP,EMCC_OP,EMCH_IP,EMCH_OP,EMCH_ST,
     &     EMDH_IP,EMDH_OP,EMDH_ST
      DATA EMCC_IP/1740.D0/,EMCC_OP/600.D0/,EMCH_IP/1161.D0/,
     &     EMCH_OP/886.D0/,EMCH_ST/3030.D0/,EMDH_IP/852.D0/,
     &     EMDH_OP/650.D0/,EMDH_ST/2224.D0/
!
!***********************************************************************
! Subroutine PAH_MODES
! Given:
!       MXMODES = dimensioning information
!       NC      = number of C atoms in PAH
!       NH      = number of H atoms (normal H) in PAH
!       ND      = number of D atoms in PAH
!
! Returns:
!       NMODES = number of vibrational degrees of freedom = 3*(NC+NH-2)
!       EMODES(J) = frequencies (cm-1) of vibrational modes J=1-NMODES
!                   in order of increasing energy
! and, via COMMON/EVCOM/,
!       EV1_COM = EMODES(1) = lowest energy mode
!
! Synthetic mode spectrum is created assuming
! NC-2 C-C out-of-plane modes distributed according to modified
!      2-dimensional Debye model with Debye frequency EMCC_OP=600 cm-1 
! 2*(NC-2) C-C in-plane modes distributed according to modified
!      2-dimensional Debye model with Debye frequency EMCC_IP=1740 cm-1
! NH C-H out-of-plane vibrations with frequency EMCH_OP=886 cm-1
! NH C-H in-plane vibrations with frequency     EMCH_IP=1161 cm-1
! NH C-H stretching modes with frequency        EMCH_ST=3030 cm-1
! ND D-H out-of-plane vibrations with frequency EMDH_OP=650 cm-1
! ND D-H in-plane vibrations with frequency     EMDH_IP=852 cm-1
! ND D-H stretching modes with frequency        EMDH_ST=2224 cm-1
!
! We use mode spectrum which is designed to resemble 2-d Debye
! mode spectrum, but adjusted so that
! * for NC.le.54 we assume a disk geometry, with lowest frequency 
!   C-C in-plane mode frequency, and lowest frequency C-C out-of-plane
!   mode frequency varying as (NC-2)**0.5
! * for 55.le.NC.le.101 we assume a transition from disk to sphere
!   geometry, with lowest frequency C-C out-of-plane and C-C in-plane
!   mode frequencies independent of NC (largest dimension of grain is
!   not changing).
! * for NC.gt.102 we assume a spherical grain geometry, with
!   lowest-frequency C-C in-plane and out-of-plane modes varying as
!   (NC-2)**(-1/3)
! We achieve this by taking (for C-C modes)
!
!                     (1-beta)
! E_k = Theta * sqrt[ --------*(k-delta_k) + beta ]
!                       N_m
!
! where
! delta_k = 0.5 for k=1 
! delta_k =  1      k=2 or 3
! delta_k = 0.5     k.ge.4
!
! Theta = Debye frequency for either C-C in-plane or C-C out-of-plane modes
!
! beta = 0 for N_C .le. 54
!
!          (N_m-N_md)
! beta = --------------   for 55 .le. N_C .le. 101
!        N_md*(2*N_m-1)
!
!            1         N_m    (102)**(2/3)
! beta = --------- * [ ---- * ------------ - 1 ]
!        (2*N_m-1)     N_md    (NC)**(2/3)
!
! where for out-of-plane modes, N_m =NC-2
!                               N_md=52
!                               N_ms=100
!           in-plane modes, N_m =2*(NC-2)
!                           N_md=104
!                           N_ms=200
!
! B.T. Draine, Princeton Univ. Obs., 00.01.26
! History:
! 99.01.26 (BTD) first written
! 99.02.09 (BTD) corrected error (C-C in-plane and C-C out-of-plane 
!                degeneracies were reversed)
! 99.02.09 (BTD) modified to use subroutine HPSORT to provide modes
!                in order of increasing frequency
! 99.02.10 (BTD) adjust distribution of lowest three modes
!                make modifications for transition from planar to
!                spherical grain.
!                We now assume planar for NC.le.52
!                              spherical for NC.ge.102
!                transition from planar to spherical for 52<NC<102
! 99.02.21 (BTD) corrected calculation of beta for NC>102
! 00.07.31 (BTD) added COMMON/EVCOM/EV1_COM to return lowest energy
!                vibrational mode
! end history
!***********************************************************************
!*** diagnostic
!      write(0,*)'pah_modes.f ckpt 0'
!      write(0,*)' MXMODES=',MXMODES
!      write(0,*)' NC=',NC
!      write(0,*)' NH=',NH
!      write(0,*)' ND=',ND
!***
      NMODES=0
      N1=NC-2
      DJOFF=0.5
!
! Compute beta for C-C out-of-plane and in-plane modes:
!
      N_MD=52
      N_S=102
      IF(NC.LE.54)THEN
         BETA=0.D0
      ELSEIF(NC.GT.54.AND.NC.LT.102)THEN
         BETA=DBLE(N1-N_MD)/DBLE(N_MD*(2*N1-1))
      ELSE
         BETA=((DBLE(N1)/N_MD)*(DBLE(N_S)/NC)**(2./3.)-1.)/DBLE(2*N1-1)
      ENDIF
!
! C-C out-of-plane modes:
!
! Adjust location of first 3 modes to bring agreement with mode spectrum
! of coronene C_{24}H_{12}
!
      EMODES(1)=SQRT((1.-BETA)*(1.-DJOFF)/DBLE(N1)+BETA)*EMCC_OP
!*** diagnostic
!      write(0,*)'pah_modes ckpt 1, emodes(1)=',emodes(1)
!***
      EMODES(2)=SQRT((1.-BETA)*(1.5-DJOFF)/DBLE(N1)+BETA)*EMCC_OP
      EMODES(3)=SQRT((1.-BETA)*(2.5-DJOFF)/DBLE(N1)+BETA)*EMCC_OP
!*** diagnostic
!      write(0,*)'pah_modes ckpt 2'
!      write(0,*)' N1=',N1
!***
      DO J=4,N1
         EMODES(J)=SQRT((1.-BETA)*(DBLE(J)-DJOFF)/DBLE(N1)+BETA)*EMCC_OP
      ENDDO
      NMODES=NMODES+N1
!
! C-C in-plane modes
!
      N1=2*(NC-2)
      EMODES(NMODES+1)=SQRT((1.-BETA)*(1.-DJOFF)/DBLE(N1)+BETA)*EMCC_IP
      EMODES(NMODES+2)=SQRT((1.-BETA)*(1.5-DJOFF)/DBLE(N1)+BETA)*EMCC_IP
      EMODES(NMODES+3)=SQRT((1.-BETA)*(2.5-DJOFF)/DBLE(N1)+BETA)*EMCC_IP
      DO J=4,N1
         EMODES(NMODES+J)=SQRT((1.-BETA)*(DBLE(J)-DJOFF)/DBLE(N1)+BETA)*
     &                    EMCC_IP
      ENDDO
      NMODES=NMODES+N1
!
! C-H modes:
!
      IF(NH.GT.0)THEN
         N1=NH
         DO J=1,N1
            EMODES(NMODES+J)=EMCH_OP
         ENDDO
         NMODES=NMODES+N1
         DO J=1,N1
            EMODES(NMODES+J)=EMCH_IP
         ENDDO
         NMODES=NMODES+N1
         DO J=1,N1
            EMODES(NMODES+J)=EMCH_ST
         ENDDO
         NMODES=NMODES+N1
      ENDIF
!
! C-D modes:
!
      IF(ND.GT.0)THEN
         N1=ND
         DO J=1,N1
            EMODES(NMODES+J)=EMDH_OP
         ENDDO
         NMODES=NMODES+N1
         DO J=1,N1
            EMODES(NMODES+J)=EMDH_IP
         ENDDO
         NMODES=NMODES+N1
         DO J=1,N1
            EMODES(NMODES+J)=EMDH_ST
         ENDDO
         NMODES=NMODES+N1
      ENDIF
!
! Now sort the modes in order of increasing frequency
!
      CALL HPSORT(NMODES,EMODES)
!
! Calculation is complete
!
      EV1_COM=EMODES(1)
      RETURN
      END










