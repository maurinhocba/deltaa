FUNCTION  vecdot(n,v1,v2)
!*****************************************************************************
!
!**** tridimensional vectorial product of two vectors  v1 x v2 -> v3
!
!*****************************************************************************
IMPLICIT NONE

  !--- Dummy variables
  REAL(kind=8) :: vecdot
  INTEGER(kind=4),INTENT(IN):: n
  REAL(kind=8),INTENT(IN) :: v1(n), v2(n)

  vecdot = DOT_PRODUCT(v1,v2)

RETURN
END FUNCTION vecdot
