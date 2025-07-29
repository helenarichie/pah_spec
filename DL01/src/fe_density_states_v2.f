      SUBROUTINE FE_DENSITY_STATES(METHOD,A,NATOM,UMAX,
     &                             NMODES,EMODES,JEMAX0,E,DE,
     &                             LNGSUM,EBS,LNG,NDIME,NDIMMODES)

      IMPLICIT NONE

!                 fe_density_states_v2

! arguments:

      CHARACTER METHOD*5
      INTEGER JEMAX0,NATOM,NDIME,NDIMMODES,NMODES
      DOUBLE PRECISION 
     &   A,DE,EBS,UMAX
      DOUBLE PRECISION
     &   E(0:NDIME),
     &   EMODES(NDIMMODES),
     &   LNG(0:NDIME),
     &   LNGSUM(0:NDIME)

! local variables:

      INTEGER I,J,JEMAX,NC,ND,NH,NM,NSET
      DOUBLE PRECISION BETA,DJOFF,EM3D,EMAX,HC,KB,XX

!-----------------------------------------------------------------------
! subroutine FE_DENSITY_STATES        v2
! given
!     NDIME     = dimensioning parameter
!     NDIMMODES = dimensioning parameter
!     METHOD    = 'stati' for fully statistical calculation
!     A         = grain radius (cm)
!     NATOM     = number of atoms in grain
!     UMAX      = max energy (erg) to consider
! returns
!     NMODES    = number of vibrational modes
!     EMODES(1-NMODES) = energies (cm-1) of vibrational modes

! and, if METHOD='stati', also return:
!     E(0-NDIME)= upper bounds of energy sub-bins
!     DE        = energy increment (cm-1) for Beyer-Swinehart algorithm
!     LNGSUM[J] = ln[gsum(JE)] where gsum(JE) = cumulative number of
!                 states up to energy E[J] = upper bound of bin J
!     EBS       = energy (cm-1) below which Beyer-Swinehart algorithm can
!                 be used without overflow
!     JEMAX0    = number of elements used in arrays E and LNGSUM

! NB: This routine should be consistent with 
!     /u/draine/work/lib/util/src/fe_modes.f
!-----------------------------------------------------------------------
! History
! 11.02.19 (BTD) * adapted from sil_density_states
! 16.03.21 (BTD) v2
! 16.03.23 (BTD) * switched to T_Debye=420K
!                simplified formula for mode distribution
!                made consistent with /u/draine/work/lib/util/src/fe_modes.f
! end history
!-----------------------------------------------------------------------

      HC=6.62607D-27*2.99792D10
      KB=1.38065D-16

! only calculate lngu[] for statistical treatments;
! in case of thermal treatments (real/Debye modes), 
! assign a number to lngu[] and return; 00.02.18.
 
      DO I=0,NDIME
         E(I)=0.D0
         LNGSUM(I)=0.D0
      ENDDO

      DO I=1,NDIMMODES
         EMODES(I)=0.D0
      ENDDO

      IF(A.LE.25.D-8)THEN

         CALL FE_MODES(NDIMMODES,NATOM,NMODES,EMODES)

      ELSEIF(A.GT.25.D-8)THEN

         NMODES=50000

! (NMODES = 50000 for a~30A; note 'stati' is limited to a<=25A)!!

! Fe Debye temperature = 420 K (Ashcroft & Mermin Table 23.3, based
! on fitting experimental specific heat to Debye model at point where
! c_v = (3/2)nk_B

         EM3D=420.D0*KB/HC
         NM=3*(NATOM-2)
         DJOFF=0.5

         DO J=1,100
            EMODES(J)=EM3D*((DBLE(J)-DJOFF)/DBLE(NM))**(1./3.)
         ENDDO

! this has now defined EMODES from J=1 to 100
! ??? why is NMODES=50000 ???

      ENDIF

      IF(METHOD.NE.'stati')THEN
         JEMAX0=0
         EBS=0.D0
         DE=1.D-2

! 08.03.23 (BTD) added following to extend calculation to larger
!                silicate grains.

         XX=EMODES(2)-EMODES(1)
         IF(XX.LT.0.02D0)THEN
            IF(XX.GT.0.01D0)THEN
               DE=5.D-3
            ELSEIF(XX.GT.0.004D0.AND.XX.LE.0.01D0)THEN
               DE=2.D-3
            ELSEIF(XX.GT.0.002D0.AND.XX.LE.0.004D0)THEN
               DE=1.D-3
            ELSEIF(XX.GT.0.001D0.AND.XX.LE.0.002D0)THEN
               DE=5.D-4
            ELSEIF(XX.GT.0.0004D0.AND.XX.LE.0.001D0)THEN
               DE=2.D-4
            ELSEIF(XX.GT.0.0002D0.AND.XX.LE.0.0004D0)THEN
               DE=1.D-4
            ELSEIF(XX.GT.0.0001D0.AND.XX.LE.0.0002D0)THEN
               DE=5.D-5
            ELSE
               WRITE(0,*)'ckpt z in fe_density_states: STOP'
               STOP
            ENDIF
         ENDIF
            
!*** diagnostic 08.03.23
!      write(0,*)'sil_dens_states ckpt 1'
!      write(0,*)'beta=',beta
!      write(0,*)'em3d=',em3d
!      write(0,*)'djoff=',djoff
!      write(0,*)'emodes(1)=',emodes(1)
!      write(0,*)'emodes(2)=',emodes(2)
!      write(0,*)'emodes(3)=',emodes(3)
!      write(0,*)'emodes(4)=',emodes(4)
!***
         RETURN
      ENDIF	

! modified on 00.07.14: let DE=1.0 or 0.1!!!!
! since nset is usually set as 20, so take 
! EMODES(20)-EMODES(19) for comparison.
! on July 19, 2000: set nset=11 (=18 modes).
! note Bruce's 10 implies 11 here!
! also note here nset=18 is actually the 18th mode
! which corresponds to the 10th bin (nset=11 in 
! sil_enthalpy_state).

! 08.03.23 (BTD) modified to keep DE small enough to avoid
!                trouble even for large grains
!                Note: stati option has not been tested with this change
      DE=1.D0
      XX=EMODES(2)-EMODES(1)
      IF(XX.LT.2.D0)THEN
         IF(XX.GT.1.D0.AND.XX.LE.2.D0)THEN
            DE=5.D-1
         ELSEIF(XX.GT.0.4D0.AND.XX.LE.1.D0)THEN
            DE=2.D-1
         ELSEIF(XX.GT.0.2D0.AND.XX.LE.0.4D0)THEN
            DE=1.D-2
         ELSEIF(XX.GT.0.1D0.AND.XX.LE.0.2D0)THEN
            DE=5.D-2
         ELSEIF(XX.GT.0.04D0.AND.XX.LE.0.1D0)THEN
            DE=2.D-2
         ELSEIF(XX.GT.0.02D0.AND.XX.LE.0.04D0)THEN
            DE=1.D-2
         ELSEIF(XX.GT.0.01D0.AND.XX.LE.0.02D0)THEN
            DE=5.D-3
         ELSEIF(XX.GT.0.004D0.AND.XX.LE.0.01D0)THEN
            DE=2.D-3
         ELSEIF(XX.GT.0.002D0.AND.XX.LE.0.004D0)THEN
            DE=1.D-3
         ELSEIF(XX.GT.0.001D0.AND.XX.LE.0.002D0)THEN
            DE=5.D-4
         ELSEIF(XX.GT.0.0004D0.AND.XX.LE.0.001D0)THEN
            DE=2.D-4
         ELSEIF(XX.GT.0.0002D0.AND.XX.LE.0.0004D0)THEN
            DE=1.D-4
         ELSEIF(XX.GT.0.0001D0.AND.XX.LE.0.0002D0)THEN
            DE=5.D-5
         ELSE
            WRITE(0,*)'ckpt z in sil_density_states: STOP'
            STOP
         ENDIF
      ENDIF

!      NSET=18
!      IF((EMODES(NSET)-EMODES(NSET-1)).GT.2.D0)THEN
!         DE=1.D0
!      ELSEIF((EMODES(NSET)-EMODES(NSET-1)).GT.1.D0 .AND.
!     &       (EMODES(NSET)-EMODES(NSET-1)).LE.2.D0)THEN
!         DE=0.5D0
!      ELSEIF((EMODES(NSET)-EMODES(NSET-1)).GT.0.4D0 .AND.
!     &       (EMODES(NSET)-EMODES(NSET-1)).LE.1.D0)THEN
!         DE=0.2D0
!      ELSEIF((EMODES(NSET)-EMODES(NSET-1)).GT.0.2D0 .AND.
!     &       (EMODES(NSET)-EMODES(NSET-1)).LE.0.4D0)THEN
!         DE=0.1D0
!      ELSEIF((EMODES(NSET)-EMODES(NSET-1)).GT.0.1D0 .AND.
!     &       (EMODES(NSET)-EMODES(NSET-1)).LE.0.2D0)THEN
!         DE=0.05D0
!      ELSEIF((EMODES(NSET)-EMODES(NSET-1)).GT.0.04D0.AND.
!     &       (EMODES(NSET)-EMODES(NSET-1)).LE.0.1D0)THEN
!         DE=0.02D0
!      ELSEIF((EMODES(NSET)-EMODES(NSET-1)).LE.0.04D0)THEN
!         WRITE(0,*)'DE needs to be refined!'
!         PAUSE
!      ENDIF

! for Fe, set NC=NH=ND=0

      NC=0
      NH=0
      ND=0

! 1: first estimate EBS;

      EMAX=-UMAX/HC

      CALL DENS_STATES(NDIME,NDIMMODES,NMODES,NC,ND,NH,NATOM,
     &                 DE,EMAX,EMODES,JEMAX,E,LNG,LNGSUM)
      EBS=EMAX*HC

! 2: then do the real job!

      EMAX=UMAX/HC
      CALL DENS_STATES(NDIME,NDIMMODES,NMODES,NC,ND,NH,NATOM,
     &                 DE,EMAX,EMODES,JEMAX,E,LNG,LNGSUM)

      JEMAX0=JEMAX+1

      RETURN
      END
