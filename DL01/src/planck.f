      FUNCTION PLANCK(WAVELENGTH,T)
      IMPLICIT NONE

c arguments:

      DOUBLE PRECISION PLANCK,T,WAVELENGTH

c local variables:

      DOUBLE PRECISION C,EKT,H,K
c-----------------------------------------------------------------------
c function PLANCK
c given
c     wavelength (cm)
c     T = temperature (K)
c returns
c                                 2*h*c^2
c     planck = B_lambda = ----------------------  (erg s-1 cm-3)
c                         lambda^5*[exp(E/kT)-1]

c              where E = h*c/lambda
c
c-----------------------------------------------------------------------
      H=6.62607D-27
      C=2.99792D10
      K=1.38065D-16
      EKT=H*C/(WAVELENGTH*K*T)
      IF(EKT.LT.100.)THEN
         PLANCK=2.*H*C*C/(WAVELENGTH**5.)/
     &          (DEXP(EKT)-1.)
      ELSE
         PLANCK=0.
      ENDIF
      RETURN
      END
