      SUBROUTINE DENS_STATES(NDIME,NDIMMODES,NMODES,NC,ND,NH,NSIL,
     &                       DE,EMAX,EMODES,JEMAX,ENERGY,LNG,LNGSUM)
      IMPLICIT NONE

! Parameters:

      INTEGER NDIME,NDIMMODES

! Arguments:

      INTEGER JEMAX,NC,ND,NH,NSIL,NMODES
      DOUBLE PRECISION DE,EMAX
      DOUBLE PRECISION 
     &   EMODES(1:NDIMMODES),ENERGY(0:NDIME),
     &   LNG(0:NDIME),LNGSUM(0:NDIME)

! Local variables:

      INTEGER DJE,J,JE,JEBS,JM
      DOUBLE PRECISION CK,DE2,DELTAE,E,EBS,EKT,FAC,TEMP

!***********************************************************************
! Subroutine DENS_STATES
!
! Given:
!
!       NDIME,NDIMMODES = dimensioning information
!       NMODES     = number of vibrational modes
!       NC         = number of C atoms
!       ND         = number of D atoms
!       NH         = number of H atoms
!       NSIL       = number of silicate atoms
!       EMODES(1-NMODES) = energies (cm-1) of vibrational modes
!       DE         = width of energy bins (cm-1)
!       EMAX       = maximum energy to consider (cm-1)
!                  < 0 to return min[E_bs,|input Emax|]
!                    where E_bs = energy (cm-1) up to which
!                    Beyer-Swinehart calculation can be done using
!                    IEEE double-precision arithmetic
! Returns:
!
! If input EMAX > 0:
!       JEMAX      = number of E>0 energy bins
!       ENERGY(JE) = energy (cm-1) of upper-limit of bin JE
!       LNG(JE)    = ln[g(JE)] where
!                    g(JE) = number of states in energy bin JE
!       LNGSUM(JE) = ln[gsum(JE)] where
!                    gsum(JE) = cumulative number of states up to and
!                               including bin JE
!                 
!                    index JE runs from 0 to JEMAX
!                    energy bin JE = [JE-1,JE]*DE
!
! If input EMAX < 0:
!    output EMAX   = min[E_bs,|input Emax|]
!                    where E_bs = energy (cm-1) up to which
!                    Beyer-Swinehart calculation can be done using
!                    IEEE double-precision arithmetic


! The Beyer-Swinehart algorithm (Beyer & Swinehart 1973; Stein & 
! Rabinovitch 1973) is used to calculate g(JE) and gsum(JE) until we 
! encounter overflow problems
! For larger energies, we integrate dS=dE/T
!
! References:
! Beyer, T., & Swinehart, D.F. 1973, Comm. of the A.C.M. 16, 379
! Stein, S.E., & Rabinovitch, B.S. 1973, J.Chem.Phys. 58, 2438-2445
!
! B.T. Draine, Princeton Univ. Observatory, 00.01.26
! History:
! 00.01.26 (BTD) written and tested using coronene
! 00.01.31 (BTD) modified to allow use for large grains by
!                using asymptotic behavior of gsum to estimate ln(gsum)
! 00.03.28 (BTD) modified to allow
!                1. both PAH and silicate grains
!                2. to be capable of computing density of states at
!                   high energies by integrating dS=dQ/T
! 00.03.31 (BTD) modified to add ENERGY array to argument list,
!                and to speed up calculation of degeneracy at high
!                energies by using logarithmically-spaced energy bins
! 00.07.12 (BTD) modified to accept EMAX < 0 as input and then to
!                determine max energy E_bs for Beyer-Swinehart
!                algorithm to be used.
! 00.07.31 (BTD) modified so that fractional bin width is not less than
!                0.01 when integrating dS=dE/T
!                (and also to indicate progress...)
! 00.10.27 (BTD) changed 0.01 to 0.01D0 in argument list of MAX()
! 00.11.21 (BTD) added sanity check
! 00.11.28 (BTD) modified to gradually increase DE after switching from
!                Beyer-Swinehart calculation to integration of dS=dQ/T
!                DE is allowed to increase by at most a factor 2 per
!                step.
!                Add term {d ln deltaE}  to calculation of {d ln g}
! 00.11.29 (BTD) modified to set JEMAX=NDIME when EMAX/DE > NDIME
! 16.03.28 (BTD) cosmetic changes
! end history
!
!***********************************************************************

      IF(EMAX.GT.0.)THEN
         JEMAX=NINT(EMAX/DE)
      ELSE

! Arrive here only when called with EMAX < 0 to test for maximum
! energy up to which Beyer-Swinehart is feasible.

         JEMAX=NINT(-EMAX/DE)
      ENDIF

! sanity check:

      IF(JEMAX.GT.NDIME)THEN
         WRITE(0,*)'Trouble in dens_states: JEMAX=',JEMAX,' > NDIME=',
     &             NDIME,' : set JEMAX to NDIME'
         JEMAX=NDIME
      ENDIF

! NB: We will use storage space LNG and LNGSUM to store g and gsum [as 
! opposed to ln(g) and ln(gsum)] while we are using the Beyer-Swinehart
! algorithm to construct g and gsum.  After we conclude with the Beyer-
! Swinehart algorithm, we will replace g and gsum by ln(g) and ln(gsum).
!
! Note that the ground state E=0 has degeneracy 1 and ln(g)=0
!
! Initialize g(J):

 0100 DO JE=1,JEMAX
         LNG(JE)=0.D0
      ENDDO
      LNG(0)=1.D0

! Execute Beyer-Swinehart algorithm:
! Note that for time being LNG is used to store g

      DO JM=1,NMODES
         DJE=NINT(EMODES(JM)/DE)

! do not allow mode frequency to be rounded to zero...

         IF(DJE.EQ.0)DJE=1
         DO JE=DJE,JEMAX

! temporarily use LNG(JE) to store G(JE)

            LNG(JE)=LNG(JE)+LNG(JE-DJE)
         ENDDO
         IF(LNG(JEMAX).GT.1.D305)THEN

! If Beyer-Swinehart algorithm is leading to overflow, try
! limiting it to a smaller maximum energy level JEMAX*DE

            IF(JM.LT.NMODES/32)THEN
               JEMAX=INT(JEMAX/2.D0)
            ELSEIF(JM.LT.NMODES/16)THEN
               JEMAX=INT(JEMAX/1.4D0)
            ELSEIF(JM.LT.NMODES/8)THEN
               JEMAX=INT(JEMAX/1.2D0)
            ELSEIF(JM.LT.NMODES/4)THEN
               JEMAX=INT(JEMAX/1.1D0)
            ELSE
               JEMAX=INT(JEMAX/1.05)
            ENDIF
            GOTO 0100
         ENDIF
      ENDDO

! Successfully completed Beyer-Swinehart calculation.

      JEBS=JEMAX
      EBS=JEMAX*DE
      IF(EMAX.LE.0.D0)THEN
         EMAX=EBS
      ENDIF
      LNGSUM(0)=LNG(0)

! LNGSUM is used temporarily to store gsum

      DO JE=1,JEMAX
         LNGSUM(JE)=LNGSUM(JE-1)+LNG(JE)
      ENDDO

! Now compute ln(g) and ln(sum(g))
! Some energy bins may have zero degeneracy.  Set degeneracy to 1e-99
! so that logarithm can be taken.

      DO JE=0,JEMAX
         ENERGY(JE)=JE*DE
         LNG(JE)=LOG(LNG(JE)+1.D-99)
         LNGSUM(JE)=LOG(LNGSUM(JE))
      ENDDO

      WRITE(0,*)'Beyer-Swinehart algorithm used up to E=',JEMAX*DE

      IF(JEMAX*DE.GE.EMAX)THEN

         RETURN
      ENDIF

! Need to extrapolate to estimate degeneracies of higher
! levels.

! We determine temperature T(E) for higher levels and
! evaluate entropy change dS=dE/T 

! Old approach:
!    use energy bins with fractional width equal to last
!    Beyer-Swinehart energy bin IF fractional width is more than 0.01
!    Otherwise set fractional width to 0.01

! New approach:
!    Use energy bins with fractional width equal to last Beyer-Swinehart
!    energy bin IF fractional width is more than 0.01
!    Otherwise increase energy bin width by factor 2 at a time
!    Until reach fractional bin width of 0.01

      DE2=DE
      JE=JEMAX
 1770 CONTINUE
      IF(DE2/ENERGY(JE).LT.0.005D0)THEN
         JE=JE+1
         DE2=2*DE2
         ENERGY(JE)=ENERGY(JE-1)+DE2
         GOTO 1770
      ENDIF

      FAC=1.01D0
! 16.04.14 (BTD) change
!      JM=NINT(LOG(EMAX/ENERGY(JE))/LOG(FAC))
      JM=INT(LOG(EMAX/ENERGY(JE))/LOG(FAC))+1

      DO J=1,JM
         ENERGY(J+JE)=ENERGY(JE)*FAC**J
      ENDDO
      JEMAX=JE+JM

! Completed construction of energy grid ENERGY(1-JM)
! Now complete calculation of ln(g) and ln(gsum)

      DO J=JEBS+1,JEMAX

! E is midpoint energy

         E=0.5*(ENERGY(J)+ENERGY(J-1))

         IF(NC.GT.0)THEN
            CALL PAH_TEMP(0,E,NC,NH,ND,TEMP,EKT,CK)
         ELSEIF(NSIL.GT.0)THEN
            CALL SIL_TEMP(0,E,NSIL,TEMP,EKT,CK)
         ENDIF
         DELTAE=ENERGY(J)-ENERGY(J-1)

! integrate dS = dE/T to calculate density of states:
! let g = degeneracy of bin, and deltaE = bin width
! then
!     S     = k ln dN/dE + const
!    dS     = k d ln (g/deltaE)
!    dE/T   = k d ln g - d ln(deltaE)
!   d ln g  = dE/kT + d ln(deltaE)

         LNG(J)=LNG(J-1)+1.48377*DELTAE/TEMP+
     &          LOG(DELTAE/(ENERGY(J-1)-ENERGY(J-2)))
         IF((LNG(J)-LNGSUM(J-1)).LT.100.)THEN
            LNGSUM(J)=LNGSUM(J-1)+LOG(1.+EXP(LNG(J)-LNGSUM(J-1)))
         ELSE
            LNGSUM(J)=LNG(J)
         ENDIF
      ENDDO

      RETURN
      END


