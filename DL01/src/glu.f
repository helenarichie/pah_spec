      SUBROUTINE GLU(WAVELENGTH,UIA,UIB,UFA,UFB,GFI)
      IMPLICIT NONE
      DOUBLE PRECISION 
     &   GFI,UFA,UFB,UIA,UIB,WAVELENGTH

c local variables

      DOUBLE PRECISION 
     &   E,E1,E2,E3,E4,DUF,DUI,HC

c ****************************************************
c - calculate Glu factor. Feb.07, 2000, Princeton.   -
c ----------------------------------------------------

c wavelength: cm; E: erg;

      HC=6.62607D-27*2.99792D10
      E=HC/WAVELENGTH

      E1=UIA-UFB
      E2=DMIN1((UIA-UFA),(UIB-UFB))
      E3=DMAX1((UIA-UFA),(UIB-UFB))
      E4=UIB-UFA

      DUI=UIB-UIA
      DUF=UFB-UFA

      IF(E1.NE.E2)THEN

         IF(E.LT.E1)THEN
            GFI=0.D0
         ELSEIF(E.GE.E1 .AND. E.LT.E2)THEN
            GFI=(E-E1)/DUF
         ELSEIF(E.GE.E2 .AND. E.LE.E3)THEN
	    IF(E2.EQ.E3)THEN
               GFI=1.
            ELSEIF(E2.NE.E3)THEN
               GFI=DMIN1(DUI,DUF)/DUF
            ENDIF
         ELSEIF(E.GT.E3 .AND. E.LT.E4)THEN
            GFI=(E4-E)/DUF
         ELSEIF(E.GE.E4)THEN
            GFI=0.D0
         ENDIF

      ELSEIF(E1.EQ.E2)THEN

         IF(E.LT.E1 .OR. E.GE.E4)THEN
            GFI=0.D0
         ELSEIF(E.GE.E1 .AND. E.LT.E4)THEN
            GFI=1.D0
         ENDIF
      ENDIF
      RETURN
      END
