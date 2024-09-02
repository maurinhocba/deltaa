SUBROUTINE invmtx(a,b,deter,nsize)
!***********************************************************************
!
!***this routine inverts a square matrix "a"
!
!     a  -  given nsize*nsize matrix
!     b  -  inverse of "a"
!     nsize  -  matrix size
!
!***********************************************************************
IMPLICIT NONE
INTEGER (kind=4),INTENT(IN) :: nsize
REAL (kind=8),INTENT(IN) :: a(nsize,nsize)
REAL (kind=8),INTENT(OUT):: b(nsize,nsize),deter
END SUBROUTINE invmtx
