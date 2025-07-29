      SUBROUTINE SIL_MODES(MXMODES,NAT,NMODES,EMODES)
      IMPLICIT NONE

! Arguments:

      INTEGER MXMODES
      INTEGER NAT,NMODES
      DOUBLE PRECISION EMODES(1:MXMODES)

! Local variables:

      INTEGER J,NM
      DOUBLE PRECISION BETA,DJOFF,EM2D,EM3D,TD2D,TD3D
      DATA TD2D/500.D0/,TD3D/1500.D0/
      SAVE TD2D,TD3D

!***********************************************************************
! Subroutine SIL_MODES
! Given:
!       MXMODES = dimensioning information
!       NAT     = number of atoms in cluster
!
! and internal variables which could be adjusted:
!       TD2D    = Debye temperature for 2-d Debye model
!       TD3D    = Debye temperature for 3-d Debye model
!
! Returns:
!       NMODES = number of vibrational degrees of freedom = 3*(NC+NH-2)
!       EMODES(J) = frequencies (cm-1) of vibrational modes J=1-NMODES
!                   in order of increasing energy
!
! Synthetic mode spectrum is created assuming
! 2*(NA-2) modes characterized by 2-d Debye density of states
!                and Debye temperature Theta=500K
! (NA-2) modes characterized by 3-d Debye density of states
!                and Debye temperature Theta=1500K
!
! We assume that lowest frequency mode for both 2-d and 3-d Debye
! mode spectrum hbar*omega_1 = k*T_Debye*N_m**(-1/3), where N_m is 
! the number of modes either the n=2 or n=3 mode spectrum.
!
! Thus hbar*omega_j = k*Theta*[(1-beta)*(j-delta_j)/N_m + beta]^{1/n}
! where n=2 or 3
! where
! delta_j = 0.5 for j=1 
! delta_j =  1      j=2 or 3
! delta_j = 0.5     j.ge.4
!
! old formula:
! For n=2, beta = [2*N_m^{1/3}-1]/(2*N_m-1)
! For n=3, beta = 1/(2*N_m-1)
!
! new formula:
! For n=2, beta = [N_m^{1/3}-1]/(2*N_m-1)
! For n=3, beta = 0
!
! B.T. Draine, Princeton Univ. Obs., 00.02.21
! History:
! 00.02.21 (BTD) adapted from pah_modes
! 00.02.22 (BTD) further modified
! 00.02.23 (BTD) complete debugging
! 00.03.05 (BTD) decided to modify formula for beta so that 
!                beta=0 for n=3
! 16.03.28 (BTD) cosmetics
! end history
!***********************************************************************
      NMODES=0
      EM2D=TD2D/1.4388D0
      EM3D=TD3D/1.4388D0

! First do 2-d Debye mode spectrum:

      NM=2*(NAT-2)
      DJOFF=0.5
      BETA=(DBLE(NM)**(1./3.)-1.D0)/(2.*NM-1.D0)

! 2/3 of modes distributed according to 2-d Debye model
! with Theta=500K

      EMODES(1)=SQRT((1.-BETA)*(1.-DJOFF)/DBLE(NM)+BETA)*EM2D
      EMODES(2)=SQRT((1.-BETA)*(2.-DJOFF-0.5)/DBLE(NM)+BETA)*EM2D
      EMODES(3)=SQRT((1.-BETA)*(3.-DJOFF-0.5)/DBLE(NM)+BETA)*EM2D
      DO J=4,NM
         EMODES(J)=EM2D*SQRT((1.-BETA)*(DBLE(J)-DJOFF)/DBLE(NM)+BETA)
      ENDDO
      NMODES=NMODES+NM

! 1/3 of modes distributed according to 3-d Debye
! model with Theta=1500K

      NM=NAT-2
      BETA=0.
      EMODES(NMODES+1)=EM3D*((1.-BETA)*(1.-DJOFF)/DBLE(NM)+
     &                       BETA)**(1./3.)
      EMODES(NMODES+2)=EM3D*((1.-BETA)*(2.-DJOFF-0.5)/DBLE(NM)+
     &                       BETA)**(1./3.)
      EMODES(NMODES+3)=EM3D*((1.-BETA)*(3.-DJOFF-0.5)/DBLE(NM)+
     &                       BETA)**(1./3.)
      DO J=4,NM
         EMODES(NMODES+J)=EM3D*
     &             ((1.-BETA)*(DBLE(J)-DJOFF)/DBLE(NM)+BETA)**(1./3.)
      ENDDO
      NMODES=NMODES+NM

! Now sort the modes in order of increasing frequency

      CALL HPSORT(NMODES,EMODES)
      RETURN
      END










