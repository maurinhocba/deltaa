SUBROUTINE vecuni(n,v,modul)
!*****************************************************************************
!
!****this routine compute de length of the vector v and
!    converts it to a unit one
!
!*****************************************************************************
IMPLICIT NONE

INTEGER (kind=4),INTENT(IN) :: n
REAL (kind=8),INTENT(IN OUT) :: v(n)
REAL (kind=8),INTENT(OUT) :: modul

modul= SQRT(DOT_PRODUCT(v,v))
IF(modul > 0)v = v/modul
RETURN

END SUBROUTINE vecuni
