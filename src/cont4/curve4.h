SUBROUTINE curve4(x,surf,np)
!
! Compute nodal normals and segment Z-increments of a surface
!
USE cont4_db
IMPLICIT NONE
INTEGER (KIND=4), INTENT (IN) :: np
REAL (KIND=8), INTENT(IN) :: x(:,:)
TYPE (surf4_db) :: surf
END SUBROUTINE curve4
