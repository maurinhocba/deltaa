SUBROUTINE inveas(n,v1,v2)
!*****************************************************************************
!
!***  vector assign:    v1(n) -> v2(n)
!
!*****************************************************************************
IMPLICIT NONE

INTEGER (kind=4),INTENT(IN) :: n,v1(n)
INTEGER (kind=4),INTENT(OUT):: v2(n)

v2 = v1
RETURN

END SUBROUTINE inveas
