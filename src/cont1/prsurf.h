      SUBROUTINE prsurf(iwric,bottom,top)

!.... write rigid surface data for visualization
!.... open file for contact forces between surfaces
!.... update flags BOTTOM and TOP

      USE cont1_db
      IMPLICIT NONE
!     dummy arguments
      INTEGER (kind=4), INTENT(IN) :: iwric
      LOGICAL, INTENT(OUT) :: bottom,top
      END SUBROUTINE prsurf
