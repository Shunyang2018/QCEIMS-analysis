      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  Program file name: maxwell.f90                                         !
      !                                                                         !
      !  © Tao Pang 2006                                                        !
      !                                                                         !
      !  Last modified: July 7, 2006                                            !
      !                                                                         !
      !  (1) This F90 program is created for the book, "An Introduction to      !
      !      Computational Physics, 2nd Edition," written by Tao Pang and       !
      !      published by Cambridge University Press on January 19, 2006.       !
      !                                                                         !
      !  (2) No warranties, express or implied, are made for this program.      !
      !                                                                         !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      MODULE CSEED
        INTEGER :: SEED
      END MODULE CSEED
      !
      ! Subroutine to draw initial velocities fron the Maxwell
      ! distribution for a given mass M, temperature T, and
      ! particle number N.
      
      SUBROUTINE MAXWELL (N, M, T, V)
      !
      !
        IMPLICIT NONE
        INTEGER, INTENT (IN) :: N
        INTEGER :: I, NV = 3*N, NG = NV-6
        REAL, INTENT (IN) :: M, T
        REAL :: V1, V2, EK, VS
        REAL, INTENT (OUT), DIMENSION (N) :: V
      !
      ! Assign a Gaussian number to each velocity component
      !
        DO I = 1, N-1, 2
          CALL RANG (V1, V2)
          V(I)   = V1
          V(I+1) = V2
        END DO
      !
      ! Scale the velocity to satisfy the partition theorem
      !
        EK = 0.0
        DO I = 1, N
          EK = EK + V(I)*V(I)
        END DO
        VS = SQRT(M*EK*NV/(NG*T))
        DO I = 1, N
          V(I) = V(I)/VS
        END DO
      END SUBROUTINE MAXWELL
      !
      SUBROUTINE RANG (X,Y)
      !
      ! Two Gaussian random numbers generated from two uniform random
      ! numbers.
      !
        IMPLICIT NONE
        REAL, INTENT (OUT) :: X,Y
        REAL :: PI,R1,R2,R,RANF
      !
        PI = 4.0*ATAN(1.0)
        R1 = -ALOG(1.0-RANF())
        R2 = 2.0*PI*RANF()
        R1 = SQRT(2.0*R1)
        X  = R1*COS(R2)
        Y  = R1*SIN(R2)
      END SUBROUTINE RANG
      !
      FUNCTION RANF() RESULT (CR)
      !
      ! Function to generate a uniform random number in [0,1]
      ! following x(i+1)=a*x(i) mod c with a=7** 5 and
      ! c=2** 31-1.  Here the seed is a global variable.
      !
        USE CSEED
        IMPLICIT NONE
        INTEGER :: H, L, T, A, C, Q, R
        DATA A/16807/, C/2147483647/, Q/127773/, R/2836/
        REAL :: CR
      !
        H = SEED/Q
        L = MOD(SEED, Q)
        T = A*L - R*H
        IF (T .GT. 0) THEN
          SEED = T
        ELSE
          SEED = C + T
        END IF
        CR = SEED/FLOAT(C)
      END FUNCTION RANF