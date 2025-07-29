      SUBROUTINE GRA_PAH_DENSITY_STATES(METHOD,A,NC,NH,ND,UMAX,
     &                                  NMODES,EMODES,JEMAX0,E,DE,
     &                                  LNGSUM,EBS,LNG,NDIME,NDIMMODES)

      IMPLICIT NONE

! arguments:

      CHARACTER METHOD*5
      INTEGER NC,NH,ND,NDIME,NDIMMODES,NMODES,JEMAX0
      DOUBLE PRECISION A,DE,EBS,UMAX
      DOUBLE PRECISION
     &   E(0:NDIME),
     &   EMODES(NDIMMODES),
     &   LNG(0:NDIME),
     &   LNGSUM(0:NDIME)

! local variables

      INTEGER I,J,JEMAX,N_MD,N_S,N1,NATOM,NSET
      DOUBLE PRECISION BETA,DJOFF,EMAX,HC

!-----------------------------------------------------------------------
! Subroutine GRA_PAH_DENSITY_STATES
!
! given:
!     METHOD = 'stati' for fully statistical calculation
!     A      = grain radius (cm)
!     NC     = number of C atoms
!     NH     = number of H atoms (not including D)
!     ND     = number of D atoms
!     UMAX   = maximum energy (erg) to consider
!     NDIME  = dimensioning for energy sub-bins E(J)
!     NDIMMODES = dimensioning for EMODES(J)

! returns:
!     NMODES = number of vibrational modes
!     EMODES(1-NMODES) = energies (cm-1) of vibrational modes

! and, if METHOD='stati', also return:
!     JEMAX0 = number of sub-bins for tabulating LNGSUM
!     E[J] = upper bound (cm-1) of sub-bin J, 1.le.J.le.JEMAX0
!     DE   = energy increment (cm-1) for Beyer-Swinehart algorithm
!     LNGSUM[J] = ln[gsum(JE)] where gsum(JE) = cumulative number of
!                 states up to energy E[J] = upper bound of bin J
!     EBS  = energy (cm-1) below which Beyer-Swinehart algorithm can
!                 be used without overflow
!

! written by Aigen Li, Princeton University
! modified by B.T. Draine, Princeton University
! uses subroutine DENS_STATES developed by B.T. Draine
! History:
! 00.11.20 (BTD) cosmetic changes
! 00.11.28 (BTD) added LNG to argument list
!                eliminate local arrays E_ORI, LNG_ORI, LNGSUM_ORI
! end history
!-----------------------------------------------------------------------
!*** diagnostic
!      write(0,*)'gra_pah_density_states ckpt 0'
!      write(0,*)'NDIMMODES=',NDIMMODES
!***
      HC=6.62607D-27*2.99792D10

! only calculate lngu[] when METHOD = 'stati'
! otherwise, assign zero to lngu[] and return
 
      DO I=0,NDIME
         E(I)=0.D0
         LNGSUM(I)=0.D0
      ENDDO

! for a<=25A, call PAH_MODES.f
! to calculate complete mode spectrum; 
! for a>25A, just
! calculate the first 100 levels.

      DO I=1,NDIMMODES
         EMODES(I)=0.D0
      ENDDO

      IF(A.LE.25.D-8)THEN

! subroutine PAH_MODES returns NMODES, EMODES
!*** diagnostic
!         write(0,*)'gra_pah_density_states ckpt 1'
!         write(0,*)' NDIMMODES=',NDIMMODES
!         write(0,*)' NC=',NC
!         write(0,*)' NH=',NH
!         write(0,*)' ND=',ND
!***
         CALL PAH_MODES(NDIMMODES,NC,NH,ND,NMODES,EMODES)
!*** diagnostic
!         write(0,*)'gra_pah_density_states ckpt 2'
!***
      ELSEIF(A.GT.25.D-8)THEN

         NMODES=50000

! (for a~30A; note 'stati' is limited to a<=25A)!!

         N1=NC-2
         DJOFF=0.5

         N_MD=52
         N_S=102
         IF(NC.LE.54)THEN
            BETA=0.D0
         ELSEIF(NC.GT.54 .AND. NC.LT.102)THEN
            BETA=DBLE(N1-N_MD)/DBLE(N_MD*(2*N1-1))
         ELSE
            BETA=((DBLE(N1)/N_MD)*(DBLE(N_S)/NC)**(2./3.)-1.)/
     &           DBLE(2*N1-1)
         ENDIF

         EMODES(1)=SQRT((1.-BETA)*(1.-DJOFF)/DBLE(N1)+BETA)*600.D0
         EMODES(2)=SQRT((1.-BETA)*(1.5-DJOFF)/DBLE(N1)+BETA)*600.D0
         EMODES(3)=SQRT((1.-BETA)*(2.5-DJOFF)/DBLE(N1)+BETA)*600.D0
         DO J=4,100
            EMODES(J)=SQRT((1.-BETA)*(DBLE(J)-DJOFF)/DBLE(N1)+BETA)*
     &                600.D0
         ENDDO

      ENDIF

! -------------------------------------------------------

      IF(METHOD.NE.'stati')THEN
         JEMAX0=0
         EBS=0.D0
         DE=1.D-2
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

      NSET=18

      IF((EMODES(NSET)-EMODES(NSET-1)).GT.2.D0)THEN
         DE=1.D0
      ELSEIF((EMODES(NSET)-EMODES(NSET-1)).GT.1.D0 .AND.
     &       (EMODES(NSET)-EMODES(NSET-1)).LE.2.D0)THEN
	 DE=0.5D0
      ELSEIF((EMODES(NSET)-EMODES(NSET-1)).GT.0.4D0 .AND.
     &       (EMODES(NSET)-EMODES(NSET-1)).LE.1.D0)THEN
         DE=0.2D0
      ELSEIF((EMODES(NSET)-EMODES(NSET-1)).GT.0.2D0 .AND.
     &       (EMODES(NSET)-EMODES(NSET-1)).LE.0.4D0)THEN
         DE=0.1D0
      ELSEIF((EMODES(NSET)-EMODES(NSET-1)).GT.0.1D0 .AND.
     &       (EMODES(NSET)-EMODES(NSET-1)).LE.0.2D0)THEN
         DE=0.05D0
      ELSEIF((EMODES(NSET)-EMODES(NSET-1)).GT.0.04D0.AND.
     &       (EMODES(NSET)-EMODES(NSET-1)).LE.0.1D0)THEN
         DE=0.02D0
      ELSEIF((EMODES(NSET)-EMODES(NSET-1)).LE.0.04D0)THEN
         WRITE(0,FMT='(A,A)')'gra_pah_density_states ckpt 3:',
     &                       'DE needs to be refined!'
         WRITE(0,FMT='(A)')' Fatal problem: Halt execution'
         STOP
      ENDIF

! natom is used for silicates.  here set to zero for graphite

      NATOM=0

! 1: first estimate EBS = maximum energy for Beyer-Swinehart algorithm
!    set EMAX to negative value so that DENS_STATES will return
!    EMAX = min[EBS,|input EMAX|]

      EMAX=-UMAX/HC
!*** diagnostic
!      write(0,*)'gra_pah_density_states ckpt 4'
!***
      CALL DENS_STATES(NDIME,NDIMMODES,NMODES,NC,ND,NH,NATOM,
     &                 DE,EMAX,EMODES,JEMAX,E,LNG,LNGSUM)
!*** diagnostic
!      write(0,*)'gra_pah_density_states ckpt 5'
!***
!-------------------------------------------
! sanity check 00.11.21

      if(jemax.gt.ndime)then
         write(0,*)'error A -- jemax > ndime in gra_pah_density_states'
         stop
      endif
!--------------------------------------------

      EBS=EMAX*HC

! 2: then do the real job!

      EMAX=UMAX/HC

!*** diagnostic
!      write(0,*)'gra_pah_density_states ckpt 6'
!***
      CALL DENS_STATES(NDIME,NDIMMODES,NMODES,NC,ND,NH,NATOM,
     &                 DE,EMAX,EMODES,JEMAX,E,LNG,LNGSUM)
!*** diagnostic
!      write(0,*)'gra_pah_density_states ckpt 7'
!***

      JEMAX0=JEMAX+1

      RETURN
      END
