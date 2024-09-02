SUBROUTINE invecz(n,v1)

!***  zero vector v1(n)

IMPLICIT NONE

INTEGER (kind=4),INTENT(IN) :: n
INTEGER (kind=4),INTENT(OUT) :: v1(n)

v1 = 0
RETURN

END SUBROUTINE invecz
