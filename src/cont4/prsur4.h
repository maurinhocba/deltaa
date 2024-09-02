      SUBROUTINE prsur4(iwric,bottom,top)

!.... write rigid surface data for visualization
!.... open file for contact forces between surfaces
!.... update flags BOTTOM and TOP

      USE cont4_db
      IMPLICIT NONE
!     dummy arguments
      INTEGER (kind=4), INTENT(IN) :: iwric
      LOGICAL, INTENT(OUT) :: bottom,top
      END SUBROUTINE prsur4
