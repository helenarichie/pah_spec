      SUBROUTINE TRAP_INTEG(N,X,Y,DELTA_LGX,Z)
      IMPLICIT NONE

! arguments

      INTEGER N
      DOUBLE PRECISION DELTA_LGX,Z
      DOUBLE PRECISION X(N),Y(N)

! local variables

      INTEGER I
!-----------------------------------------------------------------------
! subroutine TRAP_INTEG
! given
!     N = dimensions of arrays x,y
!     X[] = array
!           NOTE!! It is assumed that x is uniform in log(x)
!     Y[] = array
!     DELTA_LGX = delta(log10(x)) = increments in log10(x)

! returns
!                x[N]
!     Z   = \integral y*dx = \integral y*x*dln(x)
!                x[1]

! Written by Aigen Li, Princeton University
! History
! 00.11.13 (BTD) edited
! 16.02.19 (BTD) edited
!-----------------------------------------------------------------------
!*** diagnostic
!      write(0,*)'trap_integ ckpt 0, N=',N
!***
      Z=0.5*(X(N)*Y(N)+X(1)*Y(1))

      DO I=2,N-1
         Z=Z+X(I)*Y(I)
      ENDDO	   
      Z=Z*DLOG(10.D0)*DELTA_LGX
      RETURN
      END
