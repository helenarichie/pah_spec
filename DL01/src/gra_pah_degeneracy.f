      SUBROUTINE GRA_PAH_DEGENERACY(METHOD,A,NC,NH,ND,UMAX,NSTATE,
     &                              U,UA,UB,LNGU,NSET,JEMAX0,E,DE,
     &                              LNGSUM,EBS)

      IMPLICIT NONE

c arguments

      CHARACTER METHOD*5
      INTEGER 
     &   JEMAX0,NC,ND,NH,NSET,NSTATE
      DOUBLE PRECISION 
     &   A,DE,EBS,UMAX
      DOUBLE PRECISION
     &   E(JEMAX0),
     &   LNGSUM(JEMAX0),
     &   LNGU(NSTATE),
     &   U(NSTATE),
     &   UA(NSTATE),
     &   UB(NSTATE)

c local variables

      INTEGER I,NUA,NUB
      DOUBLE PRECISION HC,LNGSUMUA,LNGSUMUB

c*** diagnostic
      integer jj
c***
c-----------------------------------------------------------------------
c subroutine GRA_PAH_DEGENERACY
c given:
c     METHOD = 'stati'
c              (for METHOD.NE.'stati', assign LNGU[]=0 and return)
c     A      = grain radius (cm)
c     NC     = number of C atoms
c     NH     = number of H atoms (not including D)
c     ND     = number of D atoms
c     UMAX   = maximum energy (erg)
c     NSTATE = dimensioning info (max number of energy bins)
c     U[]    = midpoint bin energy (erg)
c     UA[]   = energy of lower limit of bin (erg)
c     UB[]   = energy at upper limit of bin (erg)
c     JEMAX0 = number of energy bins actually used
c     NSET   = *** not used ***
c     E[J]   = energy (cm-1)
c     LNGSUM[J] = total number of energy states with energy 
c              less than E[J] , 1.le.J.le.JEMAX0
c returns
c     LNGU[] = ln[g] , where g = degeneracy of bin 
c
c originally written by Aigen Li, Princeton Univ. Observatory
c modified by B.T. Draine
c History:
c 00.11.20 (BTD) cosmetic changes, added comments
c end history
c ------------------------------------------------------------------- 
 
      HC=6.62607D-27*2.99792D10

      IF(METHOD.NE.'stati')THEN
         DO I=1,NSTATE
            LNGU(I)=0.D0
         ENDDO
         RETURN
      ENDIF

c ground state is taken to have degeneracy = 1

      LNGU(1)=0.D0

      DO I=2,NSTATE

         IF(UA(I).LE.EBS)THEN  
            NUA=NINT(UA(I)/(HC*DE))
            LNGSUMUA=LNGSUM(NUA)
         ELSEIF(UA(I).GT.EBS)THEN
            CALL PARAB_INTERPL(E,LNGSUM,JEMAX0,UA(I)/HC,LNGSUMUA)
         ENDIF

         IF(UB(I).LE.EBS)THEN  
            NUB=NINT(UB(I)/(HC*DE))
            LNGSUMUB=LNGSUM(NUB)
         ELSEIF(UB(I).GT.EBS)THEN
            CALL PARAB_INTERPL(E,LNGSUM,JEMAX0,UB(I)/HC,LNGSUMUB)
         ENDIF

c for some cases lnGSUMua = lnGSUMub, so set lngu(i)=-100.D0
c to let gu(i)=0. (effectively)
c this allows bin to exist, but with negligible probability of
c occupation

         IF((LNGSUMUB-LNGSUMUA) .GE. 1.D-10)THEN
            LNGU(I)=LNGSUMUB+DLOG(1.-DEXP(LNGSUMUA-LNGSUMUB))
         ELSEIF((LNGSUMUB-LNGSUMUA) .LT. 1.D-10)THEN
            LNGU(I)=-100.D0

c ------------------- sanity check -----------------------------
      write(0,*)'trouble in gra_pah_degeneracy for i=',i
      write(0,*)'ub(i)=',ub(i)/HC,
     &   ' cm-1, lngsumub=',lngsumub
      write(0,*)'ua(i)=',ua(i)/HC,
     &   ' cm-1, lngsumua=',lngsumua
c --------------------------------------------------------------

            IF(LNGSUMUB.LT.LNGSUMUA)THEN
               WRITE(0,*)' in gra_pah_degeneracy, check monotonicity',
     &                   ' of lngsum from 0 to JEMAX0=',JEMAX0
               WRITE(0,*)'E(JEMAX0)=',E(JEMAX0)
               DO JJ=0,JEMAX0-1
                  IF(LNGSUM(JJ+1).LT.LNGSUM(JJ))THEN
                     WRITE(0,*)'J=',JJ,' LNGSUM(J)=',LNGSUM(JJ),
     &                         ' E=',E(JJ)
                     WRITE(0,*)'J=',JJ+1,' LNGSUM(J)=',LNGSUM(JJ+1),
     &                         ' E=',E(JJ+1)
                  ENDIF
               ENDDO
               WRITE(0,*)'------------ completed check of monotonicity',
     &                   '-- ---------'
               WRITE(0,*)' now check for strong monotonicity...'
               DO JJ=0,JEMAX0-1
                  IF(LNGSUM(JJ+1).LT.LNGSUM(JJ))THEN
                     WRITE(0,*)'J=',JJ,' LNGSUM(J)=',LNGSUM(JJ),
     &                         ' E=',E(JJ)
                     WRITE(0,*)'J=',JJ+1,' LNGSUM(J)=',LNGSUM(JJ+1),
     &                         ' E=',E(JJ+1)
                  ENDIF
               ENDDO
               WRITE(0,*)'---- completed check of strong monotonicity',
     &                   '----'

            ENDIF

c*** --------------------------------------------------------------
         ENDIF

      ENDDO

      RETURN
      END
