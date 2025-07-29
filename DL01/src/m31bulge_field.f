      SUBROUTINE M31BULGE_FIELD(WAVELENGTH,RADFLD)
      IMPLICIT NONE

c parameters

      INTEGER NTABMX
      PARAMETER(NTABMX=10000)

c arguments

      DOUBLE PRECISION WAVELENGTH,RADFLD

c common variables
c
c      INTEGER NTAB
c      DOUBLE PRECISION XT,YT
c      COMMON/TABCOM/XT(NTABMX),YT(NTABMX),NTAB

c local variables

      LOGICAL INIT
      INTEGER J,NSUM,NTAB
      DOUBLE PRECISION X,XOLD,Y,Y1,Y1SUM,Y2,Y2SUM
      DOUBLE PRECISION
     &   XT(1:NTABMX),
     &   YT(1:NTABMX)
c-----------------------------------------------------------------------
c subroutine M31BULGE_FIELD
c given
c     WAVELENGTH = wavelength (cm)
c returns
c     RADFLD = c*u_lambda (erg cm-3 s-1)
c              for the radiation field produced by the M31 bulge
c              at a distance R=1 kpc
c              as modeled by Groves et al 2012
c              arxiv 1206.2925v1 (submitted to MNRAS) Fig. 3
c              We use arithmetic mean of 
c              best fit model SED (with attenuation)
c              and unreddened stellar SED from model
c              we ignore IR emission in model by taking
c              model value at 4um and extrapolating to
c              longer wavelengths with Rayleigh-Jeans spectrum
c
c              this radiation field produces heating rate for standard
c              astrodust grain that corresponds
c              to 5.255 * MMP83 heating rate
c
c current design allows for possibility that the input table may
c have multiple entries with the same log(lambda) [as in the file
c provided by Brent Groves on 2012.10.02]
c
c B.T. Draine, Princeton University
c history
c 12.10.04 (BTD) first written
c 12.10.17 (BTD) modified to use arithmetic mean of reddened and
c                unreddened radiation fields
c end history
c
c-----------------------------------------------------------------------
      DATA INIT/.TRUE./
      SAVE INIT,NTAB,XT,YT

      IF(INIT)THEN
         OPEN(UNIT=3,FILE='/u/draine/work/irem/fields/M31_bulge.sed')

c skip  131 lines; 
c line  132 is for lambda=10.**2.9567=905 A
c line 6527 is for lambda=10.**4.6032=40087 A = 4.0087um
c for lambda > 4.0087um, we assume a Rayleigh-Jeans spectrum

         XOLD=0.
         NSUM=0
         Y1SUM=0.
         Y2SUM=0.
         NTAB=0
         DO J=1,131
            READ(3,*)
         ENDDO
         DO J=132,6527
            READ(3,*)X,Y1,Y2
            IF(X.GT.XOLD)THEN

c store previous entry

               IF(XOLD.GT.0.)THEN
                  NTAB=NTAB+1
                  XT(NTAB)=XOLD

c we take radiation field to be arithmetic mean of unreddened
c and reddened fields estimated by Groves et al

                  YT(NTAB)=LOG10(0.5*(10.**(Y1SUM/NSUM)+
     &                                10.**(Y2SUM/NSUM)))
c*** diagnostic
c                  write(0,*)'m31bulge_field_v1: ntab=',ntab,
c     &               ' xt(ntab)=',xt(ntab)
c***
               ENDIF

c start on next entry               

               XOLD=X
               Y1SUM=Y1
               Y2SUM=Y2
               NSUM=1
            ELSE
               Y1SUM=Y1SUM+Y1
               Y2SUM=Y2SUM+Y2
               NSUM=NSUM+1
            ENDIF
         ENDDO
         CLOSE(3)
         INIT=.FALSE.
      ENDIF
c-----------------------------------------------------------------------
      RADFLD=0.
      IF(WAVELENGTH.GT.9.12D-8)THEN
         X=8.+LOG10(WAVELENGTH)
         IF(X.LT.XT(NTAB))THEN

            CALL PARAB_INTERPL(XT,YT,NTAB,X,Y)

c original file contained lg(L_lambda/L_sol)*A
c for the M31 bulge stars.
c Evaluate c*u_lambda at distance R=1 kpc from this source
c L_sol/(Angstrom*4*pi*(kpc)^2)=3.198e-3 erg cm-3 s-1
c This radiation field corresponds to U=5.25

            RADFLD=3.198D-3*10.**Y
         ELSE

c assume Rayleigh-Jeans spectrum beyond longest wavelength in table
c u_lambda propto lambda^{-4}

            RADFLD=3.198D-3*10.**(YT(NTAB)-4.*(X-XT(NTAB)))
         ENDIF
      ENDIF
      RETURN
      END
