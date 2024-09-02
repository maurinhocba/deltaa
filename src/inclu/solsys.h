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
   LOGICAL, INTENT(IN) :: fact
 END SUBROUTINE Solsys
