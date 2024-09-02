SUBROUTINE proma2(a,b,c,n1,n2,n3)
!********************************************************************
!
!***this routine evaluates a matrix product
!                                                t
!           a(i,j) = b(i,k) * c(j,k)    a = b * c
!
!********************************************************************
IMPLICIT NONE

INTEGER, INTENT(IN) :: n1,n2,n3
REAL (kind=8), INTENT(IN)  :: b(n1,n3), c(n2,n3)
REAL (kind=8), INTENT(OUT) :: a(n1,n2)

a = MATMUL(b,TRANSPOSE(c))

RETURN
END SUBROUTINE proma2
