      FUNCTION provec(v1,v2)
!*****************************************************************************
!
!**** tridimensional vectorial product of two vectors  v1 x v2 -> v3
!
!*****************************************************************************
      IMPLICIT NONE

      REAL (kind=8),INTENT(IN) :: v1(3),v2(3)
      REAL (kind=8) provec(3)

      provec(1) = v1(2)*v2(3) - v1(3)*v2(2)
      provec(2) = v1(3)*v2(1) - v1(1)*v2(3)
      provec(3) = v1(1)*v2(2) - v1(2)*v2(1)
      RETURN

      END FUNCTION provec
