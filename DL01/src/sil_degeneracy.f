      SUBROUTINE SIL_DEGENERACY(METHOD,A,UMAX,NSTATE,U,UA,UB,LNGU,
     &                          NSET,JEMAX0,E,DE,LNGSUM,EBS)

      IMPLICIT NONE
      CHARACTER METHOD*5
      INTEGER JEMAX0,NSET,NSTATE

      DOUBLE PRECISION A,DE,EBS,UMAX
      DOUBLE PRECISION
     &   E(1:JEMAX0),
     &   LNGSUM(1:JEMAX0),
     &   LNGU(NSTATE),
     &   U(NSTATE),
     &   UA(NSTATE),
     &   UB(NSTATE)

      INTEGER I,NUA,NUB
      DOUBLE PRECISION HC,LNGSUMUA,LNGSUMUB

c -----------------------------------------------------------------0
c *****************************************************************
c - calculate the degeneracy of each enthalpy state of silicates. -
c - since for silicates we only consider Debye-thermal model, so  -
c - we actually do not need GU[] for silicates, thus we arbitrarily
c - set GU[]=1 for all states. March 20, 2000.                    - 
c - modified on March 24, 2000: estimate silicate degeneracies.   -
c - modified on March 29, 2000: take DENS_STATES.f out.           - 
c ----------------------------------------------------------------- 	

c lnGSUMua: exp-->{the cumulative degeneracies up to ua[]}; 
c lnGSUMub: exp-->{the cumulative degeneracies up to ub[]};
c lngu(i)=lnGSUMub[i]-lnGSUMua[i];
c EMAX: cm^{-1} (<=== umax (cm^{-1})); 
c DE: delta_E (Stein & Rabinovitch 1973) cm^{-1};
c ua(i), ub(i): erg; [/h/c] =====> cm^{-1};
c only calculate lngu[] in case of statistical treatments;
c in case of thermal treatments (real/Debye modes), assign
c a number to lngu[] and return; 00.02.18.

c History:
c 00.10.28 (BTD) modified
c 00.11.26 (BTD) cosmetic changes, added comments
c end history
c-----------------------------------------------------------------------
      HC=6.62697D-27*2.99792D10

      IF(METHOD.NE.'stati')THEN
         DO I=1,NSTATE
            LNGU(I)=0.D0
         ENDDO
         RETURN
      ENDIF
	
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
c to let gu(i)=0.! (note lnGSUMua <= lnGSUMub). 2000.02.11.

	IF((LNGSUMUB-LNGSUMUA) .GE. 1.D-10)THEN
           LNGU(I)=LNGSUMUB+DLOG(1.-DEXP(LNGSUMUA-LNGSUMUB))
	ELSEIF((LNGSUMUB-LNGSUMUA) .LT. 1.D-10)THEN
           LNGU(I)=-100.D0
	ENDIF

      ENDDO

      RETURN
      END
