      SUBROUTINE surms4(emass,x)

!.... compute nodal mass associated to surface

      USE cont4_db
      IMPLICIT NONE
!     dummy arguments
      REAL (kind=8), INTENT (IN) :: x(:,:)
      REAL (kind=8), INTENT (IN OUT) :: emass(:,:)
      END SUBROUTINE surms4
