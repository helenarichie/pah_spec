      SUBROUTINE SIMPSON_STAT_COOLING_RATE(E1,E2,E3,E4,DUI,DUF,
     &                                     NISRF,ISRF_WL,ISRF,CABS,
     &                                     EPS,RATE_COOLING)

!             simpson_stat_cooling_rate_v2

      IMPLICIT NONE

! arguments:

      INTEGER NISRF
      DOUBLE PRECISION 
     &   DUF,DUI,E1,E2,E3,E4,EPS,RATE_COOLING
      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF)

! local variables:

      INTEGER K,N
      DOUBLE PRECISION 
     &   ELOW_CUT,EUPP_CUT,H,HC,P,S1,S2,T1,T2,X

! external functions:

      DOUBLE PRECISION STAT_COOLING_FUNCT
      EXTERNAL STAT_COOLING_FUNCT

!-----------------------------------------------------------------------
! subroutine SIMPSON_STAT_COOLING_RATE
! given:
!     E1
!     E2
!     E3
!     E4
!     DUI
!     DUF
!     NISRF
!     ISRF_WL
!     ISRF
!     CABS
!     EPS

! returns:
!     RATE_COOLING

! ***************************************************************
! Calculate the cooling rate using simpson integration method -
! Originally written by Aigen Li, Princeton University
! 99.11.06 (AL)  modified
! 00.01.25 (AL)  modified
! 00.02.07 (AL)  modified
! 00.02.18 (AL)  modified
! 00.11.26 (BTD) cosmetic changes, added comments
! 18.08.13 (BTD) v2 removed A=radius from arg list(s)
! end history
!-----------------------------------------------------------------------

      HC=6.62607D-27*2.99792D10
      ELOW_CUT=HC/ISRF_WL(NISRF)
      EUPP_CUT=HC/ISRF_WL(1)

      IF(E1.LT.ELOW_CUT)E1=ELOW_CUT
      IF(E1.GT.EUPP_CUT)E1=EUPP_CUT

      IF(E2.LT.ELOW_CUT)E2=ELOW_CUT
      IF(E2.GT.EUPP_CUT)E2=EUPP_CUT

      IF(E3.LT.ELOW_CUT)E3=ELOW_CUT
      IF(E3.GT.EUPP_CUT)E3=EUPP_CUT

      IF(E4.LT.ELOW_CUT)E4=ELOW_CUT
      IF(E4.GT.EUPP_CUT)E4=EUPP_CUT

      IF(E1.EQ.E4)THEN
         RATE_COOLING=0.D0
         RETURN
      ENDIF
C --------------------------------------------------------------------  
      N=1
      H=E4-E1
      T1=0.5*H*(STAT_COOLING_FUNCT(E1,NISRF,ISRF_WL,ISRF,
     &                             CABS,E1,E2,E3,E4,DUI,DUF)+
     &          STAT_COOLING_FUNCT(E4,NISRF,ISRF_WL,ISRF,
     &                             CABS,E1,E2,E3,E4,DUI,DUF))
      S1=T1
 10   P=0.0
      DO 20 K=0,N-1
         X=E1+(K+0.5)*H
         P=P+STAT_COOLING_FUNCT(X,NISRF,ISRF_WL,ISRF,
     &                          CABS,E1,E2,E3,E4,DUI,DUF)
 20   CONTINUE
      T2=(T1+H*P)/2.0
      S2=(4.0*T2-T1)/3.0
      IF((DABS((S2-S1)/S1)).GE.EPS)THEN
         T1=T2
         N=N+N
         H=H/2.0
         S1=S2
         GOTO 10
      ENDIF
      RATE_COOLING=S2

      RETURN
      END

!***********************************************************************

      FUNCTION STAT_COOLING_FUNCT(E,NISRF,ISRF_WL,ISRF,
     &                            CABS,E1,E2,E3,E4,DUI,DUF)

      IMPLICIT NONE
      INTEGER NISRF
      DOUBLE PRECISION
     &   DUF,DUI,E,E1,E2,E3,E4,STAT_COOLING_FUNCT
      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF)

! local variables

      DOUBLE PRECISION 
     &   C,CRSSCT,GLU_E,H,HC,PI,RADFLD,RADFLD_E,WAVELENGTH

!-----------------------------------------------------------------------
! function STAT_COOLING_FUNCT
! given:
!     E
!     NISRF
!     ISRF_WL[]
!     ISRF[]
!     CABS[]
!     E1
!     E2
!     E3
!     E4
!     DUI
!     DUF

! returns

!     STAT_COOLING_FUNCT = G_lu(E)*E^3*C_abs(E)
!
! Originally written by Aigen Li, Princeton University
! History
! 00.02.07 (AL) modified
! 00.11.26 (BTD) cosmetic changes, added comments
!-----------------------------------------------------------------------

      H=6.62607D-27
      C=2.99792D10
      HC=H*C
      PI=3.14159D0

      IF(E.NE.0.)THEN

         WAVELENGTH=HC/E
         IF(WAVELENGTH.GE.ISRF_WL(1) .AND. 
     &      WAVELENGTH.LE.ISRF_WL(NISRF))THEN
            CALL PARAB_INTERPL(ISRF_WL,CABS,NISRF,WAVELENGTH,CRSSCT)
            CALL PARAB_INTERPL(ISRF_WL,ISRF,NISRF,WAVELENGTH,RADFLD)
         ELSEIF(WAVELENGTH.LT.ISRF_WL(1))THEN
            CRSSCT=CABS(1)
            RADFLD=0.
         ELSEIF(WAVELENGTH.GT.ISRF_WL(NISRF))THEN
            CRSSCT=CABS(NISRF)*(ISRF_WL(NISRF)/WAVELENGTH)**2
            RADFLD=0.
         ENDIF	

         RADFLD_E=RADFLD*WAVELENGTH**2/(HC*C)	

         IF(E1.NE.E2)THEN

            IF(E.LT.E1)THEN
               GLU_E=0.D0
            ELSEIF(E.GE.E1 .AND. E.LT.E2)THEN
               GLU_E=(E-E1)/DUF
            ELSEIF(E.GE.E2 .AND. E.LE.E3)THEN
               IF(E2.EQ.E3)THEN
                  GLU_E=1.
               ELSEIF(E2.NE.E3)THEN
                  GLU_E=DMIN1(DUI,DUF)/DUF
               ENDIF
            ELSEIF(E.GT.E3 .AND. E.LT.E4)THEN
               GLU_E=(E4-E)/DUF
            ELSEIF(E.GE.E4)THEN
               GLU_E=0.D0
            ENDIF

            STAT_COOLING_FUNCT=GLU_E*(E**3)*CRSSCT*
     &                         (8.*PI/(H*HC**2))*
     &                         (1.+HC**3/(8.*PI*E**3)*RADFLD_E)

	 ELSEIF(E1.EQ.E2)THEN

            IF(E.LT.E1 .OR. E.GE.E4)THEN
               STAT_COOLING_FUNCT=0.D0
            ELSEIF(E.GE.E1 .AND. E.LT.E4)THEN
               STAT_COOLING_FUNCT=(E**3)*CRSSCT*
     &                           (8.*PI/(H*HC**2))*
     &                           (1.+HC**3/(8.*PI*E**3)*RADFLD_E)
            ENDIF   
	 ENDIF

      ELSEIF(E.EQ.0.)THEN

	STAT_COOLING_FUNCT=0.D0

      ENDIF

      RETURN
      END

!***********************************************************************

      SUBROUTINE SIMPSON_INTRABIN_STAT_COOLING_RATE(E1,E2,DUI,NISRF,
     &                                              ISRF_WL,ISRF,CABS,
     &                                              EPS,RATE_COOLING)

      IMPLICIT NONE

! arguments

      INTEGER NISRF
      DOUBLE PRECISION DUI,E1,E2,EPS,RATE_COOLING
      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF)

! local variables

      INTEGER K,N
      DOUBLE PRECISION 
     &   ELOW_CUT,EUPP_CUT,H,HC,P,S1,S2,T1,T2,X

! external functions

      DOUBLE PRECISION INTRABIN_STAT_COOLING_FUNCT
      EXTERNAL INTRABIN_STAT_COOLING_FUNCT

!-----------------------------------------------------------------------
! subroutine SIMPSON_INTRABIN_STAT_C0OLING_RATE
! given
!     E1
!     E2
!     DUI
!     NISRF
!     ISRF_WL[]
!     ISRF[]
!     CABS[]
!     EPS

! returns
!     RATE_COOLING = power radiated (erg s-1) in intrabin transitions
!                    by grain

! We assume that the density of states is uniform within the bin.
! Then power radiated depends only on the absorption cross section
! since both initial and final micro-state will have same density
! of states.

! Originally written by Aigen Li, Princeton University Observatory
! Modified by B.T. Draine, Princeton University
! History
! 99.11.06 (AL)  First written
! 00.01.25 (AL)  Modified to avoid use of Td
! 00.02.18 (AL)  Modified to only use statistical method.
! 00.11.24 (BTD) Cosmetic changes, added comments
! 18.08.13 (BTD) v2 removed A=radius from arg list(s)
! end history
!-----------------------------------------------------------------------
      HC=6.62607D-27*2.99792D10

      ELOW_CUT=HC/ISRF_WL(NISRF)
      EUPP_CUT=HC/ISRF_WL(1)

      IF(E1.LT.ELOW_CUT)E1=ELOW_CUT
      IF(E1.GT.EUPP_CUT)E1=EUPP_CUT

      IF(E2.LT.ELOW_CUT)E2=ELOW_CUT
      IF(E2.GT.EUPP_CUT)E2=EUPP_CUT

      IF(E1.EQ.E2)THEN
         RATE_COOLING=0.D0
         RETURN
      ENDIF
C --------------------------------------------------------------------  
      N=1
      H=E2-E1
      T1=0.5*H*(INTRABIN_STAT_COOLING_FUNCT(E1,DUI,NISRF,ISRF_WL,
     &                                      ISRF,CABS)+
     &          INTRABIN_STAT_COOLING_FUNCT(E2,DUI,NISRF,ISRF_WL,
     &                                      ISRF,CABS))
      S1=T1
 10   P=0.0
      DO K=0,N-1
         X=E1+(K+0.5)*H
         P=P+INTRABIN_STAT_COOLING_FUNCT(X,DUI,NISRF,ISRF_WL,
     &                                   ISRF,CABS)
      ENDDO
      T2=(T1+H*P)/2.0
      S2=(4.0*T2-T1)/3.0
      IF((DABS((S2-S1)/S1)).GE.EPS)THEN
         T1=T2
         N=N+N
         H=H/2.0
         S1=S2
         GOTO 10
      ENDIF
      RATE_COOLING=S2

      RETURN
      END

!***********************************************************************

      FUNCTION INTRABIN_STAT_COOLING_FUNCT(E,DUI,NISRF,
     &                                     ISRF_WL,ISRF,CABS)

      IMPLICIT NONE

! arguments:

      INTEGER NISRF
      DOUBLE PRECISION 
     &   DUI,E,INTRABIN_STAT_COOLING_FUNCT
      DOUBLE PRECISION
     &   CABS(NISRF),
     &   ISRF(NISRF),
     &   ISRF_WL(NISRF)

! local variables:

      DOUBLE PRECISION 
     &   C,CRSSCT,H,HC,PI,RADFLD,RADFLD_E,WAVELENGTH

!-----------------------------------------------------------------------
! function INTRABIN_STAT_COOLING_FUNCT
! given:
!     E        = photon energy (erg)
!     DUI      = bin width (erg)
!     NISRF    = number of wavelengths for tabulated ISRF
!     ISRF_WL[]= wavelength (cm) for tabulated ISRF
!     ISRF[]   = c*u_lambda (erg cm-3 s-1) for tabulated ISRF
!     CABS[]   = C_abs (cm2) at wavelengths ISRF_WL[]

! returns

!     INTRABIN_STAT_COOLING_FUNCT = power (erg s-1 erg-1) radiated
!                    in intrabin transitions in statistical
!                    treatment of emission, with density of
!                    states approximated as constant within bin
!            8*pi         E                      (hc)^3
!        = -------* [1 - ---]*E^3*Cabs(E)*[ 1 + -------- * u_E ]
!          h^3*c^2       DUI                    8*pi*E^3

!-----------------------------------------------------------------------
! Originally written by Aigen Li, Princeton University
! Modified by B.T. Draine, Princeton University
! History
! 00.02.07 (AL)  Original version
! 00.10.30 (AL)  Statistical treatment of intrabin transitions
! 00.11.24 (BTD) Cosmetic changes, added comments
! 18.08.13 (BTD) v2 removed A=radius from arg list
! end history
!-----------------------------------------------------------------------

      H=6.62607D-27
      C=2.99792D10
      HC=H*C
      PI=3.141593D0

      IF(E.NE.0.)THEN

         WAVELENGTH=H*C/E
         IF(WAVELENGTH.GE.ISRF_WL(1) .AND. 
     &      WAVELENGTH.LE.ISRF_WL(NISRF))THEN
            CALL PARAB_INTERPL(ISRF_WL,CABS,NISRF,WAVELENGTH,CRSSCT)
            CALL PARAB_INTERPL(ISRF_WL,ISRF,NISRF,WAVELENGTH,RADFLD)
         ELSEIF(WAVELENGTH.LT.ISRF_WL(1))THEN
            CRSSCT=CABS(1)
            RADFLD=0.
         ELSEIF(WAVELENGTH.GT.ISRF_WL(NISRF))THEN
            CRSSCT=CABS(NISRF)*(ISRF_WL(NISRF)/WAVELENGTH)**2
            RADFLD=0.
         ENDIF	

         RADFLD_E=RADFLD*WAVELENGTH**2/(HC*C)	

         INTRABIN_STAT_COOLING_FUNCT=(8.*PI/(H*HC**2))*(1.-E/DUI)*
     &                               (E**3)*CRSSCT*
     &                               (1.+HC**3/(8.*PI*E**3)*RADFLD_E)

      ELSEIF(E.EQ.0.)THEN

         INTRABIN_STAT_COOLING_FUNCT=0.D0

      ENDIF
      RETURN
      END
