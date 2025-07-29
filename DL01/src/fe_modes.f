      SUBROUTINE FE_MODES(MXMODES,NAT,NMODES,EMODES)
      IMPLICIT NONE

! Arguments:

      INTEGER MXMODES
      INTEGER NAT,NMODES
      DOUBLE PRECISION EMODES(1:MXMODES)

! Local variables:

      INTEGER J,NM
      DOUBLE PRECISION BETA,DJOFF,EM2D,EM3D

!***********************************************************************
! Subroutine PAH_MODES
! Given:
!       MXMODES = dimensioning information
!       NAT     = number of atoms in cluster
!
! and internal variables which could be adjusted:
!       TD3D    = Debye temperature for 3-d Debye model
!
! Returns:
!       NMODES = number of vibrational degrees of freedom = 3*(NC+NH-2)
!       EMODES(J) = frequencies (cm-1) of vibrational modes J=1-NMODES
!                   in order of increasing energy
!
! Synthetic mode spectrum is created assuming
! 3*(NA-2) modes characterized by 3-d Debye density of states
!                and Debye temperature Theta=465K
!
! We assume that lowest frequency mode has
! hbar*omega_1 = k*T_Debye*(0.5/N_m)^(-1/3), where N_m is 
! the number of modes in the spectrum.
! and highest energy mode =
! hbar*omega_N_m = k*T_Debye*(1-0.5/N_m)^(1/3)

! B.T. Draine, Princeton Univ. Obs., 11.02.19
! History:
! 11.02.19 (BTD) adapted from sil_modes
! 16.03.23 (BTD) switched to T_Debye=420K
!                simplified formula for mode distribution
!                made consistent with fe_density_states_v2.f
! end history
!***********************************************************************
      NMODES=0
      EM3D=420./1.4388D0

! modes distributed according to 3-d Debye model with Theta=420K
! Theta=420K from Table 23.3 of Ashcroft & Mermin (1976): fitted to
! experimental specific heat at point where c_v=(3/2)nk_B

      NM=3*(NAT-2)
      DJOFF=0.5
      DO J=1,NM
         EMODES(NMODES+J)=EM3D*((DBLE(J)-DJOFF)/DBLE(NM))**(1./3.)
      ENDDO
      NMODES=NMODES+NM

! No need to sort the modes: they are already in ascending order

!      CALL HPSORT(NMODES,EMODES)

      RETURN
      END










