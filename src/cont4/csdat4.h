SUBROUTINE csdat4(maxve,iprint,label,npoin,coord,surf)

!.... read contact surface data and generate database

USE cont4_db

IMPLICIT NONE
!     arguments
INTEGER (kind=4), INTENT(IN) :: iprint,npoin,label(:),maxve
REAL (kind=8), INTENT(IN) :: coord(:,:)
TYPE (surf4_db), POINTER :: surf  !INTENT(OUT)
END SUBROUTINE csdat4
