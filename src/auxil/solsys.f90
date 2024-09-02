SUBROUTINE Solsys(A,B,n,fact)
! --------------------------------------------------
!   This SUBROUTINE solve the linear system:
!     A*x=B
!
!   The system has n equations
!   The solution x will be written in B vector
! --------------------------------------------------
  IMPLICIT NONE
  ! dummy arguments
  INTEGER(kind=4), INTENT(IN) :: n
  REAL(kind=8), INTENT(IN OUT) :: A(n,n), B(n)
  LOGICAL, INTENT (IN) :: fact
  !local variables
  INTEGER(kind=4) ::  i, j, k
  REAL(kind=8) :: L(n,n), U(n,n) !in the stack

  IF( fact )THEN
! Factorize A matrix using a LU factorization method
    DO j=1,n
      DO i=1,j-1
        DO k=i+1,n
          A(k,j)=A(k,j)-A(i,j)*A(k,i)
        END DO
      END DO
      IF ( A(j,j) == 0.0d0) THEN
        WRITE(*,*)'Termino',j,'de la Diagonal igual a cero'
        exit
       ELSE
        A(j,j)=1/A(j,j)
      END IF
      DO i=j+1,n
        A(i,j)=A(i,j)*A(j,j)
      END DO
    END DO
  END IF
! Solve the linear system
    DO i=1,n
      L(i,i)=1.0D0
      U(i,i)=A(i,i)
      DO j=i+1,n
        L(j,i)=A(j,i);  L(i,j)=0.0D0;
        U(j,i)=0.0D0;  U(i,j)=A(i,j);
      END DO
    END DO

    DO i=1,n
      DO j=1,i-1
        B(i)=B(i)-(L(i,j)*B(j))
      END DO
      B(i)=B(i)*L(i,i)
    END DO

    DO i=n,1,-1
      DO j=i+1,n
        B(i)=B(i)-(U(i,j)*B(j))
      END DO
      B(i)=B(i)*U(i,i)
    END DO

    RETURN
END SUBROUTINE Solsys
