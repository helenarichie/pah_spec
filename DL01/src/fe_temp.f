      SUBROUTINE FE_TEMP(ENERGY,NAT,TEMP,EKT,CK)

!                       fe_temp v2

      IMPLICIT NONE

! Arguments:

      INTEGER NAT
      DOUBLE PRECISION CK,EKT,ENERGY,TEMP

! Local variables:

      LOGICAL INIT
      INTEGER NMAX
      PARAMETER(NMAX=1000)
      INTEGER J,J1,JLO,JHI,NDAT
      REAL ENAT,DUM1,DUM2,X
      REAL
     &   EDAT(1:NMAX),
     &   TDAT(1:NMAX)
      DATA INIT/.true./
      SAVE EDAT,TDAT
      SAVE INIT,NDAT

!***********************************************************************
!
! Subroutine FE_TEMP
!
! Given:
!     ENERGY  = energy (cm-1)
!     NAT     = number of atoms  (must be >2)
!
! Returns:
!
!     TEMP    = temperature T (K) of Fe grain
!             = 0 if ENERGY < E_1 = first vibrationally-excited state
!     EK      = <E>/k*T , where <E>=mean thermal energy at temp. T
!     CK      = C/k  , where C(T)=d<E>/dT
!

! B.T.Draine, Princeton University Observatory
! History:
! 11.02.19 (BTD) adapted from subroutine SIL_TEMP
! 16.05.05 (BTD) v2: this is substantial change from v1
!                * eliminated argument IMODE
!                  [need to make corresponding change in subroutine
!                   fe_enthalpy_state.f]   
! end history

!***********************************************************************

      IF(INIT)THEN
         OPEN(UNIT=11,STATUS='OLD',
     &      FILE='/u/draine/work/papers/nanoFe/fe_enthalpy.dat')
         READ(11,*)
         READ(11,*)
         DO J=1,NMAX
            READ(11,*,END=1000,ERR=1000)DUM1,DUM2
            TDAT(J)=DUM1
            EDAT(J)=DUM2
         ENDDO
 1000    NDAT=J

! rezero to EDAT(1)=0:

         DO J=NDAT,1,-1
            EDAT(J)=EDAT(J)-EDAT(1)
         ENDDO

! fe_enthalpy.dat is in J/mole: convert to cm-1/atom

         DO J=1,NDAT
            EDAT(J)=(EDAT(J)*1.E7/6.022E23)/1.98645E-16
         ENDDO
         INIT=.FALSE.
      ENDIF

! first guess for array index J1

      J1=INT(NDAT/2)

! set bounds on J1

      JLO=1
      JHI=NDAT-1

! calculate internal energy per atom to compare to bulk
! use NAT-2 to omit translational and rotational degrees of freedom

      ENAT=ENERGY/REAL(NAT-2)

! check that ENAT is within tabulated range
! if not, extrapolate and return
      IF(ENAT.GE.EDAT(NDAT-1))THEN
         TEMP=TDAT(NDAT-1)+(ENAT-EDAT(NDAT-1))*
     &        (TDAT(NDAT)-TDAT(NDAT-1))/(EDAT(NDAT)-EDAT(NDAT-1))
         CK=1.4388*(EDAT(NDAT)-EDAT(NDAT-1))/(TDAT(NDAT)-TDAT(NDAT-1))
         RETURN
      ENDIF

! if within range, then do binary chop and interpolate

 2000 IF(EDAT(J1).GT.ENAT)THEN ! j1 too high
         JHI=J1
         J1=(JLO+J1)/2
      ELSE
         IF(EDAT(J1+1).LT.ENAT)THEN   ! J1 TOO LOW
            JLO=J1
            J1=(J1+JHI)/2
         ELSE                         ! j1 not too high, j1+1 not too low
            GOTO 3000
         ENDIF
      ENDIF
      GOTO 2000

 3000 X=(ENAT-EDAT(J1))/(EDAT(J1+1)-EDAT(J1))  ! linear interpolation
      TEMP=TDAT(J1)+X*(TDAT(J1+1)-TDAT(J1))    ! to find TEMP

! EKT = (energy/atom) / kT

      EKT=ENAT*1.4388/TEMP

! CK = d(energy)/dT / k

      CK=(EDAT(J1+1)-EDAT(J1))*1.4388/(TDAT(J1+1)-TDAT(J1))
      RETURN
      END
