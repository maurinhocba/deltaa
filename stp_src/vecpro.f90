SUBROUTINE vecpro(v1,v2,v3)
!*****************************************************************************
!
!**** tridimensional vectorial product of two vectors  v1 x v2 -> v3
!
!*****************************************************************************
IMPLICIT NONE

REAL (kind=8),INTENT(IN) :: v1(3),v2(3)
REAL (kind=8),INTENT(OUT):: v3(3)

v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
RETURN

END SUBROUTINE vecpro
