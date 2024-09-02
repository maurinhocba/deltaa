      SUBROUTINE surms1(emass,x)

!.... compute nodal mass associated to surface

!      USE cont1_db
      IMPLICIT NONE
!     dummy arguments
      REAL (kind=8), INTENT (IN) :: x(:,:)
      REAL (kind=8), INTENT (IN OUT)  :: emass(:,:)
      END SUBROUTINE surms1
