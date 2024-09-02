SUBROUTINE invmtx(a,b,deter,n)
!***********************************************************************
!
!***this routine inverts a square matrix "a"
!
!     a  -  given n*n matrix
!     b  -  inverse of "a"
!     n  -  matrix size
!
!***********************************************************************
IMPLICIT NONE
INTEGER (kind=4),INTENT(IN) :: n
REAL (kind=8),INTENT(IN OUT) :: a(n,n)
REAL (kind=8),INTENT(OUT):: b(n,n),deter

END SUBROUTINE invmtx
